### DESCRIPTION ################################################################
# Execute MMAA Pipeline
# Authored by; Alexander Markowitz
# Based on "DNA methylation-based classification of central nervous system tumours"
# Capper et al., nature, 2018.  PMID: 29539639

### PREAMBLE ###################################################################
library(optparse)

### PARAMETERS #################################################################
option_list <- list( 
  make_option(c("-s", "--sample-sheet"), dest = "sample_sheet", default = NULL,
              help=paste0("<string> REQUIRED: Sample sheet file: MUST BE A CSV WITH 5 COLUMNS:\n",
                          "\t\t1) sentrix barcode <string>\n",
                          "\t\t2) sample ID <string>\n",
                          "\t\t3) sample type (FFPE or Frozen) <string>\n",
                          "\t\t4) IDAT Red (Full path to red idat) <string>\n",
                          "\t\t5) IDAT Grn (Full path to green idat) <string>")),
  make_option(c("-u", "--util-dir"), default = "/usr/local/util_files/", dest = "util_dir",
              help="<string> path to util directory containing R scripts and data files [default \"%default\"]"),
  make_option(c("-C", "--CNV"), action="store_true", default=FALSE,
              help="run CNV evaluation [default = FALSE]"),
  make_option(c("-R", "--RData"), action="store_true", default=FALSE,
              help="save .RData session image at the end of analysis [default = FALSE]"),
  make_option(c("-o", "--output-dir"), default = "/docker_scratch", dest = "output_dir",
              help="<string> name of directory for html output files [default \"%default\"]")
  )

opt  <- OptionParser(option_list=option_list)
args <- parse_args(opt)

sample.file <- args$sample_sheet
util.dir    <- args$util_dir
out.dir     <- args$output_dir
run.cnv     <- args$CNV
save.image  <- args$RData

message("##############################\n",
        "SAMPLESHEET = ", sample.file, "\n",
        "UTIL DIR = ", util.dir, "\n",
        "CNV = ", run.cnv, "\n",
        "OUT DIR = ", out.dir, "\n",
        "SAVE SESSION = ", save.image, "\n",
        "##############################")

### LOAD UTILITIES #############################################################
message("Loading required packages -- ", Sys.time())
suppressMessages(suppressWarnings(suppressPackageStartupMessages(source(file.path(util.dir, "MAA_v2.R")))))
message("Done loading -- ", Sys.time())

### PERFORM CHECK SUM ##########################################################
message("Performing MD5 checksum -- ", Sys.time())
performMD5_checksum(util.dir)
message("Done loading -- ", Sys.time())

### LOAD DATA ##################################################################
source(paste0(util.dir, 'mgmtstp27/R/mgmt.R'));
load(paste0(util.dir, 'mgmtstp27/data/MGMTSTP27.rda'));
load(paste0(util.dir, 'ba_coef_v2.RData'));
load(paste0(util.dir, 'umap_object_v2.RData'));
load(paste0(util.dir, 'Mset.ctrl.RData'));
load(paste0(util.dir, 'umap_betas.RData'));
class.desc <- read.csv(paste0(util.dir, 'class_descriptions.csv'))
v12.desc   <- read.csv(paste0(util.dir, 'dkfz_ref_hierarchical_labels.csv'))
message("Done loading -- ", Sys.time())
### LOAD SAMPLE SHEET ##########################################################
sampleSheet <- read.csv(sample.file)

### RUN PIPELINE ###############################################################
message("Generating Report -- ", Sys.time())
for (i in 1:nrow(sampleSheet)) {
  
  ### SAMPLE INFO ##############################################################
  patient_info <- data.frame(
    `Sentrix Barcode` = sampleSheet$sentrix_barcode[i],
    `Patient ID` = sampleSheet$ID[i],
    `Sample Type` = sampleSheet$sample_type[i]
    )
  
  ### PROCESS IDAT FILES #######################################################
  message("Reading IDATs and classifying samples -- ", Sys.time())
  query_betas        <- suppressWarnings(idat_preprocessing(sampleSheet[i,]));
  sample_predictions <- suppressWarnings(predict_classification(sampleSheet[i,]));
  message("Done classification -- ", Sys.time())
  
  ### CNV FIT ##################################################################
  message("Running CNV analysis -- ", Sys.time())
  cnv_fit <- suppressMessages(suppressWarnings(runCNV_analysis(sampleSheet[i,])));
  message("Done CNV -- ", Sys.time())
  
  ### MGMT CLASSIFICATION ######################################################
  message("Running MGMT classification -- ", Sys.time())
  mgmt_plot <- suppressWarnings(predMGMT(sampleSheet[i,]));
  message("Done MGMT classification -- ", Sys.time())
  
  ### QC REPORT ################################################################
  message("Running QC processing -- ", Sys.time())
  qc_report <- suppressWarnings(runQC(sampleSheet[i,], query_betas));
  message("Done QC processing -- ", Sys.time())
  
  ### GENERATE REPORT ##########################################################
  message("Generating report -- ", Sys.time())
  # output
  outfile.html <- paste0(out.dir, "/", sampleSheet$ID[i], ".html")
  class_desc <- data.frame(
      class = sample_predictions[[1]]$Class,
      class_description = na.omit(class.desc$Description[class.desc$Class_v12 == sample_predictions[[1]]$Class])
  )
  message("Writing to ", outfile.html)
  rmarkdown::render(file.path(util.dir, "MAA_report.Rmd"),
                    output_file = outfile.html,
                    quiet = T)
  message("Done generating report -- ", Sys.time())
}
message("Pipeline Complete -- ", Sys.time())
### END ########################################################################