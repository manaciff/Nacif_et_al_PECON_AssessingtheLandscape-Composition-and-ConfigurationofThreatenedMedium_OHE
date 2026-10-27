# ==============================================================================
# run_all.R
# Runs the analysis chain in order
#
# Manuscript: Assessing the Landscape Composition and Configuration of
#             Threatened Medium- and Large-Sized Mammals using Habitat
#             Envelopes as Sampling Windows
# Author:     Maria Eduarda Nacif
#
# How to use
#   In the working folder of the study, open ANALISES_TCC.Rproj and press
#   Source on this file. In a clone of the repository, open R in the root of
#   the clone and run source("run_all.R"), or run Rscript run_all.R there.
#
#   Each script checks what it receives from the previous one and stops rather
#   than analyzing a reduced sample, so a failure is informative: read the
#   message before rerunning.
#
#   Scripts 00 and 01 are memory intensive. On the machine used for the study
#   they were run one at a time from a clean R session before this file, so
#   they are switched off below. Set RUN_SCRIPT_00 and RUN_SCRIPT_01 to TRUE to
#   reproduce the chain from the raw inputs in a single run.
#
#   Order matters: script 02 writes the Moran eigenvectors used by scripts 03
#   and 06, script 06 also reads the reported model of script 05 (Table 4), and
#   script 08 refits that model.
#
# Run time
#   On the machine used for the study, scripts 02 to 08 ran in about two
#   minutes. Script 05 takes the longest, because its permutation tests refit
#   the model 9,999 times. Script 00 reads 67 rasters and is slower.
# ==============================================================================

library(here)

RUN_SCRIPT_00 <- FALSE   # recompute the class metrics from the rasters
RUN_SCRIPT_01 <- FALSE   # prepare the data and run the ordinations

scripts <- c(
  if (RUN_SCRIPT_00) "00_recompute_class_metrics.R",
  if (RUN_SCRIPT_01) "01_prepare_data_and_ordination.R",
  "02_constrained_ordination_and_partitioning.R",
  "03_univariate_model_selection.R",
  "04_geometric_coupling_diagnostic.R",
  "05_ratio_artifact_test.R",
  "06_spatial_autocorrelation_checks.R",
  "08_sampling_effort_sensitivity.R"
)


# Locate this file, then the folder that holds the scripts ---------------------
# The scripts are located from the position of this file; here() is used only
# to find Dados/ and Outputs/ inside the working folder. If the position of
# this file cannot be determined (code pasted into the console), a short list
# of candidate folders is tried.
this_file_path <- function() {
  # Rscript run_all.R
  args <- commandArgs(trailingOnly = FALSE)
  hit  <- grep("^--file=", args, value = TRUE)
  if (length(hit)) return(normalizePath(sub("^--file=", "", hit[1]), mustWork = FALSE))
  # source("run_all.R") from the console or from another script
  for (i in seq_len(sys.nframe())) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(normalizePath(ofile, mustWork = FALSE))
  }
  # the Source button of RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (!is.null(p) && nzchar(p)) return(normalizePath(p, mustWork = FALSE))
  }
  NA_character_
}

self <- this_file_path()

candidates <- c(
  if (!is.na(self)) file.path(dirname(self), "R"),
  file.path(getwd(), "R"),
  file.path(getwd(), "Repositorio_GitHub", "R"),
  here("R"),
  here("Repositorio_GitHub", "R")
)
candidates <- unique(candidates[!is.na(candidates)])

has_chain <- function(p) dir.exists(p) && all(file.exists(file.path(p, scripts)))
found     <- candidates[vapply(candidates, has_chain, logical(1))]

if (!length(found)) {
  stop("The analysis scripts were not found. Folders tried:\n  ",
       paste(candidates, collapse = "\n  "),
       "\nOpen run_all.R in RStudio and press Source, or run it with ",
       "Rscript from the repository folder.", call. = FALSE)
}
script_dir <- found[1]

message("Scripts     : ", script_dir)
message("Project root: ", here())
message("Data        : ", here("Dados", "Processados"))
message("Outputs     : ", here("Outputs", "Manuscrito"))

# The repository holds a .here file, so that a clone of the repository runs on
# its own, with here() at the root of the clone. In the working folder of the
# study the repository sits one level below ANALISES_TCC.Rproj; if the R
# session starts inside the repository folder there, here() resolves to it and
# the outputs would land one level down. The check below stops only that case.
if (!file.exists(file.path(here(), "ANALISES_TCC.Rproj")) &&
    file.exists(file.path(dirname(here()), "ANALISES_TCC.Rproj"))) {
  stop("here() resolved the project root to\n  ", here(),
       "\nbut ANALISES_TCC.Rproj is in the folder above, so the chain would write ",
       "its outputs inside the repository instead of into the working folder.\n",
       "Fix: open ANALISES_TCC.Rproj and press Source on this file, or set the ",
       "working directory to the folder that holds it before sourcing.",
       call. = FALSE)
}

if (!dir.exists(here("Dados", "Processados"))) {
  stop("Folder not found: ", here("Dados", "Processados"),
       "\nhere() resolved the project root to ", here(),
       ". Open the project at the folder that holds Dados/ before running ",
       "the chain.", call. = FALSE)
}
if (!dir.exists(here("Outputs", "Manuscrito"))) {
  dir.create(here("Outputs", "Manuscrito"), recursive = TRUE)
}

# The log is written to disk line by line, so a session that stops mid-run
# still leaves a record of where it stopped.
.log_file <- file.path(here("Outputs", "Manuscrito"), "run_all_log.txt")
.say <- function(...) {
  txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", ...)
  message(txt)
  cat(txt, "\n", file = .log_file, append = TRUE, sep = "")
}
cat("", file = .log_file)
.say("run_all started")
.say("R           : ", R.version.string)
.say("working dir : ", getwd())
.say("here()      : ", here())
.say("scripts     : ", script_dir)
.say("outputs     : ", here("Outputs", "Manuscrito"))

# Modification times of every output file, to count the files each script writes.
.stamp <- function() {
  f <- list.files(here("Outputs", "Manuscrito"), recursive = TRUE, full.names = TRUE)
  stats::setNames(file.mtime(f), f)
}


# Run ----------------------------------------------------------------------------
.failed <- character(0)
for (.script in scripts) {
  message("\n", strrep("=", 78))
  .say("RUNNING  ", .script)
  message(strrep("=", 78))
  .before <- .stamp()
  .t0 <- Sys.time()
  ok <- tryCatch({
    source(file.path(script_dir, .script), echo = FALSE)
    TRUE
  }, error = function(e) {
    .say("ERROR in ", .script, ": ", conditionMessage(e))
    FALSE
  })
  .mins <- as.numeric(difftime(Sys.time(), .t0, units = "mins"))
  .after <- .stamp()
  .new_files    <- setdiff(names(.after), names(.before))
  .common_files <- intersect(names(.after), names(.before))
  .written <- length(.new_files) +
    sum(.after[.common_files] > .before[.common_files], na.rm = TRUE)
  if (ok) {
    .say(sprintf("%s finished in %.1f min, %d files written", .script, .mins, .written))
    if (.written == 0) {
      .say("WARNING: ", .script, " wrote no file. Check here() and the output path.")
    }
  } else {
    .failed <- c(.failed, .script)
    .say(sprintf("%s stopped after %.1f min", .script, .mins))
  }
  .say(sprintf("memory in use: %.0f MB", sum(gc()[, 2])))
}
if (length(.failed)) {
  .say("chain incomplete. Scripts that stopped: ", paste(.failed, collapse = ", "))
} else {
  .say("all scripts finished")
}

# The sample size is read from the file the chain uses, not asserted here.
raw <- file.path(here("Dados", "Processados"), "Data_Raw_WithCoords.csv")
n_units <- if (file.exists(raw)) nrow(read.csv(raw)) else NA_integer_
.say(sprintf("Sampling units in Data_Raw_WithCoords.csv: %s",
             ifelse(is.na(n_units), "file not found", n_units)))
message("The reported model (Table 4) is in PartV/TableS_PartV_Uncoupled_Model.csv.")
message("The full a priori model (Table 3) is in Table_S_PartIII_ConfirmatoryModel.csv.")
message("docs/TABLE_MAP.md maps every table and figure to the file that holds it.")
