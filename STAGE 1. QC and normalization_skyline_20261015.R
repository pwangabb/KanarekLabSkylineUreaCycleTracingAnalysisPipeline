# ==============================================================================
# MASTER PIPELINE: Data Import -> Split -> Non-Tracer Cleanup -> Tracer Cleanup
# Features:
#    1. Auto-Increment Numbering (01, 02, 03...)
#    2. Interactive File Selection
#    3. Smart Filtering based on Experiment Type
#    4. Double Tracer Logic (Formatting & Database Check)
#    5. Internal Standard Filter: Removes 13C, 15N, and D (Deuterium)
#    6. Nitrogen-Only Filter: Removes metabolites with ONLY N0 rows in N-Tracer mode.
#
# CHANGES IN THIS VERSION (see "FIX:" comments below):
#    A. Reference-database lookup for the no-nitrogen check (Double Tracer mode)
#       now reports every metabolite it could NOT find in the database, instead
#       of silently skipping it. Previously, a name mismatch meant the ".N0"
#       suffix was never stripped and IsoCorrectoR would abort with a syntax
#       error ("LabN0 doesn't match the required syntax") with no indication
#       of which metabolite or why.
#    C. The reference database can now be selected as either .csv or .xlsx.
#       The original script only ever used read.csv() on it, which silently
#       produces garbage ("embedded nulls" / "invalid input") if you pick an
#       actual Excel (.xlsx) file, since that's a binary format, not text.
#       File type is now detected from the extension and read accordingly
#       (requires the 'readxl' package for .xlsx).
#
#    NA-to-zero handling is UNCHANGED from the original script, per request.
# ==============================================================================

library(dplyr)
library(stringr)
library(tibble)
library(readxl)   # FIX C: needed to read .xlsx reference databases

# Helper Function: Get Next Step Number for file naming
get_next_id <- function(d) {
  files <- list.files(d, pattern="^[0-9]+")
  nums <- as.numeric(gsub("^([0-9]+)_.*", "\\1", files))
  nums <- nums[!is.na(nums)]
  if(length(nums)==0) return("01")
  return(sprintf("%02d", max(nums)+1))
}

# ==============================================================================
# PART 1: Setup, Input & Split
# ==============================================================================
cat("\n=======================================================\n")
cat("PART 1: Data Import & Configuration\n")
cat("=======================================================\n")

# 1. Select File
cat("Please select your input file (e.g., SamplesInColumns-Daily.csv)...\n")
raw_file_path <- file.choose()
source_directory <- dirname(raw_file_path)
setwd(source_directory)
cat("Working Directory set to:", source_directory, "\n")

# 2. Create Output Folder
output_folder_name <- "Research Output"
if (!dir.exists(output_folder_name)) dir.create(output_folder_name)
base_dir <- file.path(source_directory, output_folder_name)

# 3. Read Data
cat("Reading data...\n")
raw_data <- read.csv(raw_file_path, stringsAsFactors = FALSE, na.strings = c("#N/A", "NA", ""), check.names = FALSE)
colnames(raw_data) <- gsub(" Area", "", colnames(raw_data))

# 4. Define Helper Function
select_cols_by_rules <- function(all_cols, user_input) {
  selected_cols <- c("Molecule")
  rules <- trimws(strsplit(user_input, ",")[[1]])
  rules <- unlist(strsplit(rules, " "))
  rules <- rules[rules != ""]

  for (rule in rules) {
    if (grepl("-", rule)) {
      range_vals <- as.numeric(unlist(strsplit(rule, "-")))
      for (col in all_cols) {
        num_match <- regmatches(col, regexpr("\\d+", col))
        if (length(num_match) > 0 && as.numeric(num_match) >= range_vals[1] && as.numeric(num_match) <= range_vals[2]) {
          selected_cols <- c(selected_cols, col)
        }
      }
    } else {
      matches <- grep(rule, all_cols, value = TRUE)
      selected_cols <- c(selected_cols, matches)
    }
  }
  return(unique(selected_cols))
}

# 5. User Inputs
cat("\n--- Defining Sample Groups ---\n")
input_tracer <- readline(prompt = "Input for Tracer Samples (e.g., 108-143): ")
cat("\nIMPORTANT: Please include Pooled Samples in the Non-Tracer group for QC.\n")
input_nontracer <- readline(prompt = "Input for Non-Tracer Samples (e.g., 200-250 pool): ")

cat("\n--- Defining Experiment Type ---\n")
cat("1. Carbon Tracer (C)\n")
cat("2. Nitrogen Tracer (N)\n")
cat("3. Double Tracer (C + N)\n")
exp_choice <- readline(prompt = "Enter number (1-3): ")

# 6. Set Regex Patterns
m0_match_pattern <- switch(exp_choice,
                           "1" = "_C0$",
                           "2" = "_N0$",
                           "3" = "(_C0[._]?N0)$"
)

suffix_remove_pattern <- switch(exp_choice,
                                "1" = "_C[0-9]+$",
                                "2" = "_N[0-9]+$",
                                "3" = "_C[0-9]+[._]?N[0-9]+$"
)

if (is.null(m0_match_pattern)) stop("Invalid experiment choice. Please run the script again.")

# 7. Split & Save (Generating ID 01)
data_tracer <- raw_data[, select_cols_by_rules(colnames(raw_data), input_tracer), drop = FALSE]
data_nontracer <- raw_data[, select_cols_by_rules(colnames(raw_data), input_nontracer), drop = FALSE]

prefix1 <- get_next_id(base_dir) # 01
write.csv(data_tracer, file.path(base_dir, paste0(prefix1, "_Input_Tracer.csv")), row.names = FALSE)
write.csv(data_nontracer, file.path(base_dir, paste0(prefix1, "_Input_NonTracer.csv")), row.names = FALSE)

cat(sprintf("Step 1 Complete. Raw split files saved with prefix [%s].\n", prefix1))


# ==============================================================================
# PART 2: Non-Tracer Cleanup (Unchanged)
# ==============================================================================
cat("\n=======================================================\n")
cat("PART 2: Processing Non-Tracer Data\n")
cat("=======================================================\n")

df_nt <- data_nontracer
df_nt[is.na(df_nt)] <- 0

df_nt_filtered <- df_nt %>% filter(grepl(m0_match_pattern, Molecule))
df_nt_renamed <- df_nt_filtered %>%
  mutate(Molecule = str_remove(Molecule, suffix_remove_pattern))

num_cols_nt <- colnames(df_nt)[-1]
df_nt_clean <- df_nt_renamed %>%
  group_by(Molecule) %>%
  summarise(across(all_of(num_cols_nt), sum, na.rm = TRUE)) %>%
  ungroup()

row_sums_nt <- rowSums(df_nt_clean[, -1])
df_nt_final <- df_nt_clean[row_sums_nt > 0, ]

df_for_t <- as.data.frame(df_nt_final)
rownames(df_for_t) <- df_for_t$Molecule
df_for_t <- df_for_t[, -1]
df_nt_trans <- as.data.frame(t(df_for_t))
df_nt_trans <- rownames_to_column(df_nt_trans, "Sample")

prefix2 <- get_next_id(base_dir) # 02
file_nt_out <- paste0(prefix2, "_Data_NonTracer_Cleaned_Transposed.csv")
write.csv(df_nt_trans, file.path(base_dir, file_nt_out), row.names = FALSE)
cat(sprintf("Step 2 Complete. Saved: %s\n", file_nt_out))


# ==============================================================================
# PART 3: Tracer Cleanup (Added N-Tracer Filtering)
# ==============================================================================
cat("\n=======================================================\n")
cat("PART 3: Processing Tracer Data\n")
cat("=======================================================\n")

df_tr <- data_tracer
metabo_col <- colnames(df_tr)[1]
df_tr[is.na(df_tr)] <- 0

# 1. Remove Internal Standards
is_std <- grepl("13C|15N|D[0-9]+", df_tr[[1]])
if (sum(is_std) > 0) {
  df_tr <- df_tr[!is_std, ]
  cat(sprintf("Removed %d internal standard rows.\n", sum(is_std)))
}

# 2. Smart M+0 Check & Cleanup
df_check <- df_tr %>%
  mutate(BaseName = str_remove(!!sym(metabo_col), suffix_remove_pattern))

df_check$IsM0 <- grepl(m0_match_pattern, df_check[[metabo_col]])
sample_cols_tr <- colnames(df_tr)[-1]
df_check$RowSum <- rowSums(df_check[, sample_cols_tr])

# [LOGIC FOR NITROGEN TRACER]
if (exp_choice == "2") {
  cat("Filtering: Removing metabolites that only contain N0 (no 15N enrichment detected)...\n")

  labeled_mets <- df_check %>%
    filter(!grepl("_N0$", !!sym(metabo_col)) & RowSum > 0) %>%
    pull(BaseName) %>%
    unique()

  valid_mets <- df_check %>%
    filter(IsM0 == TRUE & RowSum > 0 & BaseName %in% labeled_mets) %>%
    pull(BaseName) %>%
    unique()
} else {
  valid_mets <- df_check %>%
    filter(IsM0 == TRUE & RowSum > 0) %>%
    pull(BaseName) %>%
    unique()
}

cat(sprintf("Final count: %d valid metabolites kept in tracer file.\n", length(valid_mets)))

df_tr_final <- df_check %>%
  filter(BaseName %in% valid_mets) %>%
  select(-BaseName, -IsM0, -RowSum)

# 3. Double Tracer Processing (If applicable)
if (exp_choice == "3") {
  cat("\nFormatting names to _C#.N# style...\n")
  df_tr_final[[metabo_col]] <- gsub("(_C[0-9]+)_(N[0-9]+)", "\\1.\\2", df_tr_final[[metabo_col]])

  cat("!!! IMPORTANT ACTION REQUIRED: SELECT SKYLINE DATABASE (.CSV or .XLSX) !!!\n")
  ref_file_path <- file.choose()

  # FIX C: detect file type from extension and read accordingly, instead of
  # always using read.csv(). Previously, selecting an .xlsx file here passed
  # the raw binary file straight to read.csv(), which does not understand
  # the Excel format at all -- it would just read garbage bytes as if they
  # were text, producing exactly the kind of errors you saw ("embedded
  # nulls", "invalid input found on input connection", and the downstream
  # "missing value where TRUE/FALSE needed" error when R tried to make sense
  # of the resulting nonsense column names).
  ref_ext <- tolower(tools::file_ext(ref_file_path))

  if (ref_ext %in% c("xlsx", "xls")) {
    ref_df <- as.data.frame(read_excel(ref_file_path, sheet = 1), stringsAsFactors = FALSE)
    colnames(ref_df) <- as.character(colnames(ref_df))
  } else {
    ref_df <- tryCatch({
      read.csv(ref_file_path, stringsAsFactors = FALSE, check.names = FALSE)
    }, warning = function(w) {
      read.csv(ref_file_path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-16LE")
    }, error = function(e) { stop("Could not read file.") })
  }

  formula_col <- grep("Precursor Formula", colnames(ref_df), value = TRUE, ignore.case = TRUE)[1]
  name_col <- grep("Molecule Name|Peptide|Precursor Name", colnames(ref_df), value = TRUE, ignore.case = TRUE)[1]

  if (is.na(formula_col) || is.na(name_col)) {
    stop(sprintf(
      "Could not find the expected columns in the reference database.\nLooking for a 'Precursor Formula' column and a 'Molecule Name'/'Peptide'/'Precursor Name' column.\nColumns found: %s",
      paste(colnames(ref_df), collapse = ", ")
    ))
  }

  # FIX A: trim whitespace on both sides before matching, so invisible
  # leading/trailing spaces (common after Excel round-trips) don't cause a
  # false "not found".
  ref_df[[name_col]] <- trimws(ref_df[[name_col]])

  df_tr_final$TempBaseName <- sub("_C[0-9]+\\.N[0-9]+$", "", df_tr_final[[metabo_col]])
  unique_mets <- unique(trimws(df_tr_final$TempBaseName))
  mets_without_N <- c()
  mets_not_found <- c()   # FIX A: collect every metabolite the lookup could not find

  for (met in unique_mets) {
    match_row <- ref_df[ref_df[[name_col]] == met, ]
    if (nrow(match_row) > 0) {
      if (!grepl("N(?![a-z])", match_row[[formula_col]][1], perl = TRUE)) {
        mets_without_N <- c(mets_without_N, met)
      }
    } else {
      # FIX A: previously this case was silently ignored, so the ".N0"
      # suffix was never stripped and IsoCorrectoR aborted later with an
      # opaque "LabN0 doesn't match the required syntax" error with no clue
      # as to which metabolite or why. Now we collect and report these.
      mets_not_found <- c(mets_not_found, met)
    }
  }

  if (length(mets_not_found) > 0) {
    cat(sprintf(
      "\n!!! WARNING: %d metabolite(s) were NOT found in the reference database by exact name match.\n",
      length(mets_not_found)
    ))
    cat("    Their '.N0' suffix was NOT stripped, even if they truly contain no nitrogen.\n")
    cat("    If any of these are nitrogen-free compounds, this WILL cause IsoCorrectoR to\n")
    cat("    abort with a 'LabN0 doesn't match the required syntax' error. Check these names\n")
    cat("    against the database for typos, renames, or formatting differences:\n")
    for (met in sort(unique(mets_not_found))) cat("      -", met, "\n")
    cat("\n")
  }

  df_tr_final <- df_tr_final %>%
    mutate(!!sym(metabo_col) := ifelse(TempBaseName %in% mets_without_N, sub("\\.N0$", "", !!sym(metabo_col)), !!sym(metabo_col))) %>%
    select(-TempBaseName)
}

# 4. Save Final Output
row_sums_tr <- rowSums(df_tr_final[, sample_cols_tr])
df_tr_final <- df_tr_final[row_sums_tr > 0, ]
colnames(df_tr_final)[1] <- "Measurements/Samples"

prefix3 <- get_next_id(base_dir) # 03
file_tr_out <- paste0(prefix3, "_Data_Tracer_Cleaned.csv")
write.csv(df_tr_final, file.path(base_dir, file_tr_out), row.names = FALSE)
cat(sprintf("Step 3 Complete. Saved: %s\n", file_tr_out))

cat("\n=======================================================\n")
cat("ALL DONE! Files saved in: ", base_dir, "\n")
