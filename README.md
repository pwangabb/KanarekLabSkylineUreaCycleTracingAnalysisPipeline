# Skyline Urea Cycle, NEAA, and PPP Stable Isotope Tracing Analysis Pipeline

This repository contains the full pipeline for processing ¹³C and ¹⁵N stable isotope tracing data in urea cycle metabolism, non-essential amino acid (NEAA) metabolism, and the pentose phosphate pathway (PPP), from Skyline peak picking through QC, normalization, and natural abundance correction. Three tracing experiments are supported:

1. L-Citrulline (ureido-¹³C) tracing of the urea cycle (¹³C)
2. Ammonium acetate (¹⁵N) tracing of the urea cycle and NEAAs (¹⁵N)
3. [1,2-¹³C₂]Glucose tracing of the PPP (¹³C)

---

## Pipeline Overview

```
Raw LC-MS data
      |
      v
Stage 0: Peak picking in Skyline (Database 1, 2, or 3)
      |
      v
Stage 1: QC and normalization in R (STAGE 1. QC and normalization_skyline_20261015.R)
      |
      v
Stage 2: Natural abundance correction with IsoCorrectoR (MoleculeInfo + ElementInfo files)
      |
      v
Corrected isotopologue distributions (fractions) for urea cycle metabolites, NEAAs, or PPP metabolites
```

---

## Repository Contents

| File | Description |
|------|-------------|
| `Database 1_L-Citrulline (ureido-¹³C)_urea cycle_tracing.csv` | Skyline transition list for ¹³C urea cycle tracing: citrulline (C0–C6) and argininosuccinate (C0–C10), 18 transitions |
| `Database 2_ammonium acetate (¹⁵N)_urea cycle_NEAA_tracing.csv` | Skyline transition list for ¹⁵N urea cycle and NEAA tracing: 12 metabolites, 35 transitions (all N isotopologues per metabolite; aspartate in both polarities) |
| `Database 3_Glucose 12-13C2_PPP_tracing.xlsx` | Skyline transition list for [1,2-¹³C₂]glucose PPP tracing: C0–C2 isotopologues for glycolytic and PPP intermediates, 42 transitions |
| `MoleculeInfo 1_IsoCorrectoR_L-Citrulline (ureido-¹³C)_urea cycle.csv` | IsoCorrectoR molecule file for experiment 1 |
| `MoleculeInfo 2_IsoCorrectoR_ammonium acetate (¹⁵N)_urea cycle_NEAA.csv` | IsoCorrectoR molecule file for experiment 2 |
| `MoleculeInfo 3_IsoCorrectoR_Glucose (1,2-13C2)_PPP.csv` | IsoCorrectoR molecule file for experiment 3 |
| `ElementInfo 1_IsoCorretoR_L-Citrulline (ureido-¹³C)_urea cycle.csv` | IsoCorrectoR element file for experiment 1 (¹³C tracer purity 99%) |
| `ElementInfo 2_IsoCorretoR_ammonium acetate (¹⁵N)_urea cycle_NEAA.csv` | IsoCorrectoR element file for experiment 2 (¹⁵N tracer purity 98%) |
| `ElementInfo 3_IsoCorretoR_Glucose (1,2-13C2)_PPP.csv` | IsoCorrectoR element file for experiment 3 (¹³C tracer purity 99%) |
| `STAGE 1. QC and normalization_skyline_20261015.R` | R script for data import, sample splitting, internal standard removal, and tracer data cleanup (Stage 1) |
| `STAGE 1. QC and normalization_skyline_20261015.docx` | User guide for the Stage 1 script |
| `STAGE 2_Natural abundance correction using El Maven and IsoCorrectoR.docx` | Step-by-step guide for Stage 2 |

Use the Database, MoleculeInfo, and ElementInfo files that share the same number (1, 2, or 3) for a given experiment.

---

## Stage 0: Peak Picking in Skyline

**Input:** Raw LC-MS data files
**Reference file:** The Skyline database for your experiment (`Database 1`, `Database 2`, or `Database 3`)

1. Import the transition list (`Insert > Transition List`) for your experiment into a new Skyline document.
2. Import your raw LC-MS data files.
3. Check and curate peak integration for every metabolite and isotopologue. Explicit retention times and windows are provided in each database.
4. Export a consolidated `.csv` of peak areas (Molecule names in column 1, one column per sample or injection).

**Isotopologue naming:** Databases 1 and 3 use `_C0`, `_C1`, ... (¹³C), Database 2 uses `_N0`, `_N1`, ... (¹⁵N). These suffixes are used by the Stage 1 script to recognize unlabeled (M+0) rows and to strip isotopologue labels, so do not rename them. For double-tracer experiments, name isotopologues `_C#_N#` in Skyline (for example `glutamine_C0_N1`); Stage 1 converts these to the `_C#.N#` format that IsoCorrectoR requires.

---

## Stage 1: QC and Normalization

**Script:** `STAGE 1. QC and normalization_skyline_20261015.R`

### Prerequisites
- R and RStudio
- Packages: `dplyr`, `stringr`, `tibble`, `readxl`

```r
install.packages(c("dplyr", "stringr", "tibble", "readxl"))
```

### Running the script
1. Open the script in RStudio and run it (Cmd + Option + R on Mac, Ctrl + Alt + R on Windows).
2. Select your consolidated Skyline `.csv` when prompted. The working directory is set to that file's folder, and a `Research Output` folder is created.
3. Define sample groups by entering sample numbers or keywords (for example `108-143`):
   - **Tracer samples**
   - **Non-tracer samples** (include pooled samples here for QC)
4. Select the experiment type: `1` Carbon tracer, `2` Nitrogen tracer, or `3` Double tracer (C + N).

The script splits the data into tracer and non-tracer groups, removes internal standards (rows containing ¹³C, ¹⁵N, or D labels), and keeps only metabolites that have signal. In nitrogen-tracer mode, metabolites that only have N0 signal (no ¹⁵N enrichment detected) are also removed. In double-tracer mode you are prompted to select the Skyline database (`.csv` or `.xlsx`), which is used to remove the `.N0` suffix from nitrogen-free metabolites. Any metabolite not found in the database is reported in the console.

### Output files
Files are numbered automatically and saved in `Research Output/`.

| File | Description |
|------|-------------|
| `01_Input_Tracer.csv` | Raw data for tracer samples |
| `01_Input_NonTracer.csv` | Raw data for non-tracer samples (including pools) |
| `02_Data_NonTracer_Cleaned_Transposed.csv` | Cleaned M+0 data for non-tracer samples, with samples as rows |
| `03_Data_Tracer_Cleaned.csv` | **Final tracer data** with the IsoCorrectoR measurement file header (`Measurements/Samples` in cell A1); blank out zeros before Stage 2 |

---

## Stage 2: Natural Abundance Correction

Full details are in `STAGE 2_Natural abundance correction using El Maven and IsoCorrectoR.docx`.

### Measurement file
- Row 1, Column A must read `Measurements/Samples`.
- Sample names must be unique.
- Column A lists each metabolite and isotopologue, for example `citrulline_C0`, `citrulline_C1` (¹³C) or `glutamate_N0`, `glutamate_N1` (¹⁵N). For double tracers, separate the two labels with a period (`glutamine_C0.N1`).
- Replace `0` values with blank cells.
- Save as `.csv`. `03_Data_Tracer_Cleaned.csv` from Stage 1 already has the correct header and naming, but zeros are not blanked by the script, so replace them before running IsoCorrectoR.

### Molecule file
Use the `MoleculeInfo` file for your experiment. Column A is the molecule name (case sensitive, and it must match the measurement file), column B is the formula with the tracer element and maximum number of labeled atoms (for example `C6H13N3O3LabC6` for citrulline with ¹³C, `C5H10N2O3LabN2` for glutamine with ¹⁵N), and column C is left blank for MS1 data. Add new metabolites to the file before running a new dataset.

### Element file
Use the `ElementInfo` file for your experiment. It lists natural isotope abundances, the tracer isotope mass shift, and tracer purity (¹³C 99% for experiments 1 and 3, ¹⁵N 98% for experiment 2). Update the purity value if your tracer lot differs.

### Running IsoCorrectoR
Install IsoCorrectoR from Bioconductor once (see the IsoCorrectoR documentation for the GUI package), then open **R (not RStudio)** and run:

```r
library(IsoCorrectoRGUI)
IsoCorrectionGUI()
```

Provide the paths to your measurement, molecule, and element files and an output directory, then check these settings:

- ✓ **High Resolution Mode** (required for double-tracer data)
- ✓ **Correct Tracer Impurity**
- ✓ **Calculate Mean Enrichment**

The corrected isotopologue fractions are written to the output directory.

---

## Dependencies Summary

| Tool | Used for |
|------|----------|
| Skyline | Stage 0 peak picking |
| R / RStudio | Stage 1 script |
| `dplyr`, `stringr`, `tibble`, `readxl` | Stage 1 R packages |
| R (console) | Running IsoCorrectoRGUI |
| `IsoCorrectoR`, `IsoCorrectoRGUI` | Stage 2 natural abundance correction |

## References

- IsoCorrectoR documentation: http://bioconductor.org/packages/devel/bioc/vignettes/IsoCorrectoR/inst/doc/IsoCorrectoR.html
