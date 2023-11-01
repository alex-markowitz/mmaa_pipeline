### DESCRIPTION ################################################################
# Mutli-modal Methylation Array Analysis Pipeline (MMAA)

## PREAMBLE ####################################################################
suppressMessages(library(minfi));
suppressMessages(library(randomForest));
suppressMessages(library(glmnet));
suppressMessages(library(plotly));
suppressMessages(library(conumee));
suppressMessages(library(knitr));
suppressMessages(library(tools));
suppressMessages(library(rmarkdown));
suppressMessages(library(umap));
suppressMessages(library(plotly));
suppressMessages(library(htmltools));
suppressMessages(library(rjson));
suppressMessages(library(caret));
suppressMessages(library(IlluminaHumanMethylationEPICv2manifest));
suppressMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38));
suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19));
suppressMessages(library(limma));
suppressMessages(library(optparse));
suppressMessages(library(multimode));
suppressMessages(library(EnvStats));
options(scipen = 999);

### PREPROCESSING FXNS #########################################################
readidat <- function(sampleSheet) {
    rgSet <- read.metharray(paste0('/docker_scratch/', sampleSheet$idat_grn));
    if (rgSet@annotation[1] == "Unknown") {
        rgSet@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = '20a1.hg38');
        }
    return(rgSet)
}            

MNPnormalize.illumina.control <- function (rgSet, ref=10000) {
    ## This function returns an rgset, not a methylset
    ## code duplication
    Green <- getGreen(rgSet)
    Red   <- getRed(rgSet)    

    if (minfi:::.is450k(rgSet) || minfi:::.isEPIC(rgSet) || rgSet@annotation[1] == "IlluminaHumanMethylationEPICv2") {
        AT.controls <- getControlAddress(rgSet, controlType = c("NORM_A", "NORM_T"))
        CG.controls <- getControlAddress(rgSet, controlType = c("NORM_C", "NORM_G"))
    }
    if (minfi:::.is27k(rgSet)) {
        AT.controls <- getControlAddress(rgSet, controlType = "Normalization-Red")
        CG.controls <- getControlAddress(rgSet, controlType = "Normalization-Green")
    }
    Green.avg    <- colMeans(Green[CG.controls, , drop = FALSE])
    Red.avg      <- colMeans(Red[AT.controls, , drop = FALSE])
    Green.factor <- ref/Green.avg
    Red.factor   <- ref/Red.avg
    Green        <- sweep(Green, 2, FUN = "*", Green.factor)
    Red          <- sweep(Red, 2, FUN = "*", Red.factor)
    assay(rgSet, "Green") <- Green
    assay(rgSet, "Red")   <- Red
    rgSet
}

MNPpreprocessIllumina <- function (rgSet, bg.correct = TRUE, normalize = c("controls", "no"), ref = 10000) {
    minfi:::.isRGOrStop(rgSet)
    normalize <- match.arg(normalize)

    if (normalize == "controls") {
        rgSet <- MNPnormalize.illumina.control(rgSet,ref)
    }

    if (bg.correct) {
        
        minfi:::.isRGOrStop(rgSet)
        Green <- getGreen(rgSet)
        Red   <- getRed(rgSet)
        if (minfi:::.is450k(rgSet) || minfi:::.isEPIC(rgSet) || rgSet@annotation[1] == "IlluminaHumanMethylationEPICv2") {
            NegControls <- getControlAddress(rgSet, controlType = "NEGATIVE")
        }
        if (minfi:::.is27k(rgSet)) {
            NegControls <- getControlAddress(rgSet, controlType = "Negative")
        }
        Green.bg <- apply(Green[NegControls, , drop = FALSE], 2, function(xx) {
            sort(as.vector(xx))[31]
        })
        Red.bg <- apply(Red[NegControls, , drop = FALSE], 2, function(xx) {
            sort(as.vector(xx))[31]
        })
        Green <- pmax(sweep(Green, 2, Green.bg), 0)
        Red   <- pmax(sweep(Red, 2, Red.bg), 0)
        assay(rgSet, "Green") <- Green
        assay(rgSet, "Red")   <- Red
        rgSet
    }

    out <- preprocessRaw(rgSet)
    preprocess <- sprintf("Illumina_mnp, bg.correct = %s, normalize = %s, refIntensity = %d",
                          bg.correct, normalize, ref)

    ## TODO:
    ## The manifest package version is currently not updated since 'packageVersion(getManifest(rgSet))' fails.          
    ## packageVersion expects a string
    out@preprocessMethod <- c(rg.norm = preprocess,
                              minfi = as.character(packageVersion("minfi")),
                              manifest = as.character(packageVersion(minfi:::.getManifestString(rgSet@annotation))))
    #packageVersion(getManifest(rgSet)))
    out
}
### PREPROCESSING IDAT FILES ###################################################
idat_preprocessing <- function(sampleSheet) {

    maav2.probes     <- read.csv(paste0(util.dir, "maav2.probes.csv"));
    rgSet            <- readidat(sampleSheet);
    Mset             <- MNPpreprocessIllumina(rgSet);
    row.names(Mset)  <- gsub("_.*", "", row.names(Mset));
    Mset             <- Mset[maav2.probes[,1] ,];

    # Extract methylation and unmethylation data
    query_methy   <- getMeth(Mset);
    query_unmethy <- getUnmeth(Mset);

    if (sampleSheet$sample_type == 'FFPE') {
      # Adjust methylation values of the query sample
      adjusted_methy <- log2(query_methy + 1) + methy.coef[["FFPE"]];

      # Adjust unmethylation values of the query sample
      adjusted_unmethy <- log2(query_unmethy + 1) + unmethy.coef[["FFPE"]];
    }

    else if (sampleSheet$sample_type == 'KRYO' | sampleSheet$sample_type == 'Frozen') {
      # Adjust methylation values of the query sample
      adjusted_methy <- log2(query_methy + 1) + methy.coef[["Frozen"]];

      # Adjust unmethylation values of the query sample
      adjusted_unmethy <- log2(query_unmethy + 1) + unmethy.coef[["Frozen"]];
    }

    # Apply any necessary thresholding or constraints
    adjusted_methy[adjusted_methy < 0]     <- 0;
    adjusted_unmethy[adjusted_unmethy < 0] <- 0;

    # Convert back to the original scale (if required)
    adjusted_methy   <- 2^adjusted_methy;
    adjusted_unmethy <- 2^adjusted_unmethy;

    # Recalculate query_betas, illumina-like
    query_betas <- adjusted_methy / (adjusted_methy + adjusted_unmethy + 100);
    query_betas <- as.data.frame(t(query_betas));

    return(query_betas);
}
### RUN RF CLASSIFIER ##########################################################
predict_classification <- function(sampleSheet) {

    v12.desc <- read.csv(paste0(util.dir, 'dkfz_ref_hierarchical_labels.csv'));

    for (i in 1:nrow(sampleSheet)) {

        rgSet              <- readidat(sampleSheet[i,]);
        superfamily_probes <- read.csv(paste0(util.dir, 'superfamily_cpg_probes.csv'));
        Mset               <- MNPpreprocessIllumina(rgSet);
        row.names(Mset)    <- gsub("_.*", "", row.names(Mset));
        Mset               <- Mset[superfamily_probes[,1],];

        # Extract methylation and unmethylation data
        query_methy   <- getMeth(Mset);
        query_unmethy <- getUnmeth(Mset);

        load(paste0(util.dir, "maav2_superfamily_model_data.RData"));

        if (sampleSheet$sample_type[i] == 'FFPE') {
            # Adjust methylation values of the query sample
            adjusted_methy   <- log2(query_methy + 1) + methy.coef[["FFPE_SF"]];

            # Adjust unmethylation values of the query sample
            adjusted_unmethy <- log2(query_unmethy + 1) + unmethy.coef[["FFPE_SF"]];
        }

        else if (sampleSheet$sample_type[i] == 'KRYO' | sampleSheet$sample_type[i] == 'Frozen') {
            # Adjust methylation values of the query sample
            adjusted_methy <- log2(query_methy + 1) + methy.coef[["Frozen_SF"]];

            # Adjust unmethylation values of the query sample
            adjusted_unmethy <- log2(query_unmethy + 1) + unmethy.coef[["Frozen_SF"]];
        }

        # Apply any necessary thresholding or constraints
        adjusted_methy[adjusted_methy < 0]     <- 0;
        adjusted_unmethy[adjusted_unmethy < 0] <- 0;

        # Convert back to the original scale (if required)
        adjusted_methy   <- 2^adjusted_methy;
        adjusted_unmethy <- 2^adjusted_unmethy;

        # Recalculate query_betas, illumina-like
        query_betas <- adjusted_methy / (adjusted_methy + adjusted_unmethy + 100);
        query_betas <- as.data.frame(t(query_betas));

        # Super Family Predictions and Score Calibration
        superfamily.scores         <- predict(superfamily_rf_model, query_betas, type = "prob")
        superfamily.cal.scores     <- predict(superfamily.cv.calfit, superfamily.scores, type = "response", s = superfamily.cv.calfit$lambda.1se)
        knn.superfamily.scores     <- predict(knn.superfamily.v12, query_betas, type = 'prob')
        superfamily.knn.cal.scores <- predict(knn.superfamily.v12.cv, knn.superfamily.scores, type = "response", s = knn.superfamily.v12.cv$lambda.1se)

        rf.superfamily.v12.max  <- apply(superfamily.cal.scores, 1, function(x) max(x))
        rf.superfamily.v12      <- apply(superfamily.cal.scores, 1, function(x) colnames(superfamily.cal.scores)[which.max(x)])

        knn.superfamily.v12.max      <- apply(superfamily.knn.cal.scores, 1, function(x) max(x))
        knn.superfamily.v12.name     <- apply(superfamily.knn.cal.scores, 1, function(x) colnames(superfamily.knn.cal.scores)[which.max(x)])

        query_superfamily <- ifelse(rf.superfamily.v12.max >= knn.superfamily.v12.max, rf.superfamily.v12, knn.superfamily.v12.name)

        if (length(unique(v12.desc$Family_code[v12.desc$Superfamily_code == query_superfamily])) == 1) {
            rf.family.v12.max  <- rf.superfamily.v12.max
            rf.family.v12      <- v12.desc$Family_code[v12.desc$Superfamily_code == query_superfamily][1]

            knn.family.v12.max     <- knn.superfamily.v12.max
            knn.family.v12.name    <- v12.desc$Family_code[v12.desc$Superfamily_code == query_superfamily][1]

            rf.class.v12.max       <- rf.superfamily.v12.max
            knn.class.v12.max      <- knn.superfamily.v12.max
            rf.class.v12           <- v12.desc$Class_v12_code[v12.desc$Superfamily_code == query_superfamily][1]
            knn.class.v12.name     <- v12.desc$Class_v12_code[v12.desc$Superfamily_code == query_superfamily][1]

            query_family <- ifelse(rf.family.v12.max >= knn.family.v12.max, rf.family.v12, knn.family.v12.name)

        }

        else {
            load(paste0(util.dir, "maav2_predict_family_model_data_", query_superfamily,".RData"))

            rgSet <- readidat(sampleSheet[i,]);
            family_probes <- read.csv(paste0(util.dir, query_superfamily ,"_cpg_probes.csv"));
            Mset  <- MNPpreprocessIllumina(rgSet);
            row.names(Mset) <- gsub("_.*", "", row.names(Mset));
            Mset <- Mset[family_probes[,1],];

            # Extract methylation and unmethylation data
            query_methy   <- getMeth(Mset);
            query_unmethy <- getUnmeth(Mset);

            if (sampleSheet$sample_type[i] == 'FFPE') {
                # Adjust methylation values of the query sample
                adjusted_methy <- log2(query_methy + 1) + methy.coef[["FFPE"]];

                # Adjust unmethylation values of the query sample
                adjusted_unmethy <- log2(query_unmethy + 1) + unmethy.coef[["FFPE"]];
            }

            else if (sampleSheet$sample_type[i] == 'KRYO' | sampleSheet$sample_type[i] == 'Frozen') {
                # Adjust methylation values of the query sample
                adjusted_methy <- log2(query_methy + 1) + methy.coef[["Frozen"]];

                # Adjust unmethylation values of the query sample
                adjusted_unmethy <- log2(query_unmethy + 1) + unmethy.coef[["Frozen"]];
            }

            # Apply any necessary thresholding or constraints
            adjusted_methy[adjusted_methy < 0]     <- 0;
            adjusted_unmethy[adjusted_unmethy < 0] <- 0;

            # Convert back to the original scale (if required)
            adjusted_methy   <- 2^adjusted_methy;
            adjusted_unmethy <- 2^adjusted_unmethy;

            # Recalculate query_betas, illumina-like
            query_betas <- adjusted_methy / (adjusted_methy + adjusted_unmethy + 100);
            query_betas <- as.data.frame(t(query_betas));

            # Family Predictions and Score Calibration
            family.scores          <- predict(family_rf_model, query_betas, type = "prob");
            family.cal.scores      <- predict(family.cv.calfit, family.scores, type = "response", s = family.cv.calfit$lambda.1se);
            knn.family.scores      <- predict(knn.family.v12, query_betas, type = 'prob');
            family.knn.cal.scores  <- predict(knn.family.v12.cv, knn.family.scores, type = "response", s = knn.family.v12.cv$lambda.1se);

            rf.family.v12.max  <- apply(family.cal.scores, 1, function(x) max(x));
            rf.family.v12      <- apply(family.cal.scores, 1, function(x) colnames(family.cal.scores)[which.max(x)]);

            knn.family.v12.max   <- apply(family.knn.cal.scores, 1, function(x) max(x));
            knn.family.v12.name  <- apply(family.knn.cal.scores, 1, function(x) colnames(family.knn.cal.scores)[which.max(x)]);

            query_family <- ifelse(rf.family.v12.max >= knn.family.v12.max, rf.family.v12, knn.family.v12.name);
        }

        if (length(unique(v12.desc$Class_v12_code[v12.desc$Family_code == query_family])) == 1) {
            rf.class.v12.max       <- rf.family.v12.max;
            knn.class.v12.max      <- knn.family.v12.max;
            rf.class.v12           <- v12.desc$Class_v12_code[v12.desc$Family_code == query_family][1];
            knn.class.v12.name     <- v12.desc$Class_v12_code[v12.desc$Family_code == query_family][1];
        }
        
        else {
            load(paste0(util.dir, "maav2_predict_class_model_data_", query_family,".RData"));

            rgSet           <- readidat(sampleSheet[i,]);
            class_probes    <- read.csv(paste0(util.dir, query_family ,"_cpg_probes.csv"));
            Mset            <- MNPpreprocessIllumina(rgSet);
            row.names(Mset) <- gsub("_.*", "", row.names(Mset));
            Mset            <- Mset[class_probes[,1],];

            # Extract methylation and unmethylation data
            query_methy   <- getMeth(Mset);
            query_unmethy <- getUnmeth(Mset);
            
            if (sampleSheet$sample_type[i] == 'FFPE') {
                # Adjust methylation values of the query sample
                adjusted_methy <- log2(query_methy + 1) + methy.coef[["FFPE"]];
                
                # Adjust unmethylation values of the query sample
                adjusted_unmethy <- log2(query_unmethy + 1) + unmethy.coef[["FFPE"]];
                }
            
            else if (sampleSheet$sample_type[i] == 'KRYO' | sampleSheet$sample_type[i] == 'Frozen') {
                # Adjust methylation values of the query sample
                adjusted_methy <- log2(query_methy + 1) + methy.coef[["Frozen"]];
                
                # Adjust unmethylation values of the query sample
                adjusted_unmethy <- log2(query_unmethy + 1) + unmethy.coef[["Frozen"]];
                }

            # Apply any necessary thresholding or constraints
            adjusted_methy[adjusted_methy < 0]     <- 0;
            adjusted_unmethy[adjusted_unmethy < 0] <- 0;

            # Convert back to the original scale (if required)
            adjusted_methy   <- 2^adjusted_methy;
            adjusted_unmethy <- 2^adjusted_unmethy;

            # Recalculate query_betas, illumina-like
            query_betas <- adjusted_methy / (adjusted_methy + adjusted_unmethy + 100);
            query_betas <- as.data.frame(t(query_betas));

            # Class Predictions and Score Calibration
            class.v12.scores       <- predict(class_rf_model, query_betas, type = "prob");
            class.v12.cal.scores   <- predict(class.cv.calfit, class.v12.scores, type = "response", s = class.cv.calfit$lambda.1se);
            knn.class.v12.scores   <- predict(knn.class.v12, query_betas, type = 'prob');
            class.knn.cal.scores   <- predict(knn.class.v12.cv, knn.class.v12.scores, type = "response", s = knn.class.v12.cv$lambda.1se);

            rf.class.v12.max  <- apply(class.v12.cal.scores, 1, function(x) max(x));
            knn.class.v12.max <- apply(class.knn.cal.scores, 1, function(x) max(x));
            rf.class.v12      <- apply(class.v12.cal.scores, 1, function(x) colnames(class.v12.cal.scores)[which.max(x)]);
            knn.class.v12.name     <- apply(class.knn.cal.scores, 1, function(x) colnames(class.knn.cal.scores)[which.max(x)]);

            query_class <- ifelse(rf.class.v12.max >= knn.class.v12.max, rf.class.v12, knn.class.v12.name);
            }

        calibrated_results <- data.frame(
            Sample_Name = colnames(rgSet),
            superfamily = ifelse(rf.superfamily.v12.max >= knn.superfamily.v12.max, rf.superfamily.v12, knn.superfamily.v12.name),
            superfamily.score = ifelse(rf.superfamily.v12.max >= knn.superfamily.v12.max, rf.superfamily.v12.max, knn.superfamily.v12.max),
            family = ifelse(rf.family.v12.max >= knn.family.v12.max, rf.family.v12, knn.family.v12.name),
            family.score = ifelse(rf.family.v12.max >= knn.family.v12.max, rf.family.v12.max, knn.family.v12.max),
            class.v12 = ifelse(rf.class.v12.max >= knn.class.v12.max, rf.class.v12, knn.class.v12.name),
            class.score.v12 = ifelse(rf.class.v12.max >= knn.class.v12.max, rf.class.v12.max, knn.class.v12.max)
            );

        v12_info <- data.frame(
            `Superfamily` = v12.desc$Superfamily[v12.desc$Superfamily_code == calibrated_results$superfamily[1]][1],
            `Superfamily Score` = calibrated_results$superfamily.score,
            `Family` = v12.desc$Family[v12.desc$Family_code == calibrated_results$family[1]][1],
            `Family Score` = calibrated_results$family.score,
            `Class` = v12.desc$Class_v12[v12.desc$Class_v12_code == calibrated_results$class.v12[1]][1],
            `Class Score` = calibrated_results$class.score.v12
            );
        }
    return(list(v12_info));
  }
### RUN UMAP ANALYSIS ##########################################################
generate_umap <- function(query_betas) {
    # Perform UMAP on training betas data
    embedding_train_df <- data.frame(UMAP1 = umap_result$layout[, 1], UMAP2 = umap_result$layout[, 2], Sample_Class = betas$Class);

    # Perform UMAP on query sample data
    embedding_query    <- predict(umap_result, query_betas);
    embedding_query_df <- data.frame(UMAP1 = embedding_query[, 1], UMAP2 = embedding_query[, 2]);

    # Create UMAP scatter plot
    fig <- plot_ly();

    # Get unique labels and generate a color map
    labels <- unique(embedding_train_df$Sample_Class);
    colors <- rainbow(length(labels));

    # Add training samples to the plot
    for (i in seq_along(labels)) {
        label <- labels[i]
        color <- colors[i]
        mask <- embedding_train_df$Sample_Class == label
        fig <- add_trace(fig,
                         x = ~UMAP1,
                         y = ~UMAP2,
                         data = subset(embedding_train_df, mask),
                         type = "scatter",
                         mode = "markers",
                         name = label,
                         marker = list(size = 6, color = color)
        );
    }

    # Add query sample to the plot with a text label and line pointer
    fig <- add_trace(fig,
                     x = ~UMAP1,
                     y = ~UMAP2,
                     data = embedding_query_df,
                     type = "scatter",
                     mode = "markers+text",
                     name = "Query",
                     text = "Query Sample",
                     textposition = "top center",
                     hoverinfo = "skip",
                     marker = list(size = 10, color = "red")
        );

    # Update plot layout
    fig <- layout(fig,
                  title = "UMAP Clustering",
                  xaxis = list(title = "UMAP1"),
                  yaxis = list(title = "UMAP2"),
                  showlegend = FALSE,
                  legend = list(orientation = "h", yanchor = "bottom", y = 1.02, xanchor = "right", x = 1)
        );
    return(fig);
}

### RUN MGMT CLASSIFIER ########################################################
predMGMT <- function(sampleSheet) {
    rgSet <- readidat(sampleSheet);
    MGMT_Mset            <- preprocessRaw(rgSet);
    row.names(MGMT_Mset) <- gsub("_.*", "", row.names(MGMT_Mset));
    MGMT_Mset            <- MGMT_Mset[c("cg12981137", "cg12434587"),];
    MGMT_Mval            <- log2((getMeth(MGMT_Mset)+1)/(getUnmeth(MGMT_Mset)+1));
    MGMT_Mval            <- as.data.frame(t(MGMT_Mval));
    pred                 <- MGMTpredict(MGMT_Mval);

    cutoff <- 0.358 # as per PMID 22810491
    plot <- ggplot(pred, aes(x = pred, xmin = lower, xmax = upper, y = 0)) +
        geom_point(size = 3) +
        geom_errorbar(width = 0.25) +
        xlim(c(0, 1)) +
        ylim(c(-0.75, 0.75)) +
        xlab("MGMT Score") +
        geom_vline(xintercept = cutoff, linetype = "dashed", color = "red") +
        theme_bw() +
        theme(axis.ticks.y = element_blank(),
              axis.title.y = element_blank(), 
              axis.text.y = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.x = element_line(linetype = "dotted", color = "gray65"),
              axis.text.x = element_text(face = "bold", size = 15),
              axis.title.x = element_text(face = "bold", size = 15));
    return(plot);
}

### RUN CNV ANALYSIS ###########################################################
runCNV_analysis <- function(sampleSheet) {
    rgSet              <- readidat(sampleSheet);
    Mset               <- MNPpreprocessIllumina(rgSet);

    if (Mset@annotation[1] == "IlluminaHumanMethylationEPICv2") {
        row.names(Mset)    <- gsub("_.*", "", row.names(Mset));
    }

    Mset               <- Mset[match(row.names(Mset.ctrl), row.names(Mset)),];
    colData(Mset.ctrl) <- NULL;
    Mset$type          <- "tumor";
    Mset.ctrl$type     <- "normal";
    Mset.combined      <- cbind(Mset, Mset.ctrl);

    data(exclude_regions);
    data(detail_regions);

    anno <- CNV.create_anno(array_type = "450k", 
                            exclude_regions = exclude_regions, 
                            detail_regions = detail_regions);

    anno@probes    <- anno@probes[names(anno@probes) %in% rownames(Mset.combined)];
    cnv.data       <- CNV.load(Mset.combined);
    minfi.controls <- pData(Mset.combined)$type == "normal";

    # run CNV
    cnv.fits <- list();

    for (i in 1:ncol(Mset)) {

        fit     <- CNV.fit(query = cnv.data[names(cnv.data)[!minfi.controls]][i], ref = cnv.data[names(cnv.data)[minfi.controls]], anno);
        bin     <- CNV.bin(fit);
        detail  <- CNV.detail(bin);
        segment <- CNV.segment(detail);
        segment@bin[["ratio"]][is.na(segment@bin[["ratio"]])] <- 0;

        sample.list <- list(
            fit,
            bin,
            detail,
            segment
        );
        names(sample.list)   <- c('fit', 'bin', 'detail', 'seg');
        cnv.fits[[i]]        <- sample.list;
        names(cnv.fits)[[i]] <- sampleSheet$ID[i]
    }
    return(cnv.fits);
}

### RUN QC ANALYSIS ############################################################
runQC <- function(sampleSheet, query_betas) {
    rgSet  <- readidat(sampleSheet);
    Mset   <- MNPpreprocessIllumina(rgSet);
    Rset   <- ratioConvert(Mset, what = "both", keepCN = TRUE);
    GRset  <- mapToGenome(Rset);
    rm(Rset);

    #  Calculate median intensity of X and Y chrom.
    CN     <- getCN(GRset);
    xIndex <- which(seqnames(GRset) == "chrX");
    yIndex <- which(seqnames(GRset) == "chrY");
    rm(GRset);

    # Sex Chrm Signal Med Intensity
    xMed <- colMedians(CN, rows = xIndex, na.rm = TRUE);
    yMed <- colMedians(CN, rows = yIndex, na.rm = TRUE);

    methyl   <- getMeth(Mset);
    unmethyl <- getUnmeth(Mset);

    median.df <- data.frame(
        methyl_log2_median = apply(methyl, 2, function(x) median(log2(x))),
        unmethyl_log2_median = apply(unmethyl, 2, function(x) median(log2(x)))
        );

    median.df$medIntensityOutlierDetected <- ifelse(
        median.df$methyl_log2_median < 11 |
            median.df$unmethyl_log2_median < 11,
        TRUE,
        FALSE
        );

    pop.results <- NULL;
    query_betas <- as.data.frame(t(query_betas))

    for (i in c(1:ncol(query_betas))) {

        temp.results         <- locmodes(query_betas[,i], mod0 = 2, display = FALSE, lowsup = 0, uppsup = 1);
        temp.lower.beta.mode <- temp.results$locations[1];
        temp.upper.beta.mode <- temp.results$locations[3];
        temp.anti.mode       <- temp.results$locations[2];

        temp.lower.density <- temp.results$fvalue[1];
        temp.upper.density <- temp.results$fvalue[3];
        temp.anti.density  <- temp.results$fvalue[2];

        temp.critical.bandwidth <- temp.results$cbw$bw;

        sample.results <- data.frame(
            sample.id       = colnames(query_betas)[i],
            lower.beta.mode = temp.lower.beta.mode,
            upper.beta.mode = temp.upper.beta.mode,
            anti.mode       = temp.anti.mode,
            lower.density   = temp.lower.density,
            upper.density   = temp.upper.density,
            crit.bw         = temp.critical.bandwidth
        );

        pop.results <- rbind(pop.results, sample.results);
    }

    # calculate the detection p-values
    detP           <-  detectionP(rgSet);
    detPercentage  <- 1 - colMeans(detP);
    detThres       <- colMeans(detP) > 0.05;

    # QC Report
    qc <- data.frame(
        detP = detPercentage,
        detThres = detThres,
        beta.thes = pop.results$upper.beta.mode < 0.8,
        methyl_log2_median = apply(methyl, 2, function(x) median(log2(x))) < 11,
        unmethyl_log2_median = apply(unmethyl, 2, function(x) median(log2(x))) < 11
      );
    return(qc);
}