## code to prepare `mlnhanes` dataset goes here
# ============================================================
# NHANES (2007–2018) downloader + cleaner + outcomes
# - Outcome 1: Diabetes (HbA1c / fasting glucose / self-report / meds)
# - Outcome 2: Hypertension (measured BP / self-report / meds)
#
# Author: (you)
# Requires internet for download
# ============================================================

# ---- Packages ----

library(nhanesA)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(readr)

# ---- Cycles ----
# NHANES file suffix letters:
# 2007–2008: _E, 2009–2010: _F, 2011–2012: _G, 2013–2014: _H, 2015–2016: _I, 2017–2018: _J
cycle_suffix <- c("E","F","G","H","I","J")
cycle_label  <- c("2007-2008","2009-2010","2011-2012","2013-2014","2015-2016","2017-2018")
names(cycle_suffix) <- cycle_label

# Base tables you want every cycle (best-effort)
base_stems <- c(
  "DEMO", "BMX", "BPX",
  "GHB", "GLU",
  "INS", "VID",
  "TRIGLY", "TCHOL", "HDL", "BIOPRO",
  "DIQ", "BPQ", "SMQ", "ALQ", "PAQ",
  "SLQ"
)

# CRP changed across cycles:
# - early cycles: CRP_E / CRP_F exist :contentReference[oaicite:4]{index=4}
# - later cycles: HSCRP_I / HSCRP_J exist :contentReference[oaicite:5]{index=5}
crp_stem_by_suffix <- list(
  E = "CRP",
  F = "CRP",
  G = "CRP",    # may or may not exist; safe_nhanes() will skip if absent
  H = "CRP",    # same
  I = "HSCRP",
  J = "HSCRP"
)

# Build the actual table names per cycle, e.g., "DEMO_E"
tables_for_cycle <- function(suf) {
  stems <- base_stems
  # add the correct CRP/HSCRP stem for that cycle
  stems <- c(stems, crp_stem_by_suffix[[suf]])
  unique(paste0(stems, "_", suf))
}

# ---- What to download (fairly large but still manageable) ----
# Core:
# DEMO: demographics + weights
# BMX: body measures
# BPX: blood pressure exam
#
# Labs (widely used, mostly stable naming):
# GHB: HbA1c (LBXGH)
# GLU: fasting glucose subsample (LBXGLU)
# TRIGLY: triglycerides
# TCHOL: total cholesterol
# HDL: HDL cholesterol
# CRP: C-reactive protein
# BIOPRO: standard biochemistry profile (includes creatinine, ALT/AST, etc. depending on cycle)
#
# Questionnaires:
# DIQ: diabetes questionnaire (self-report + insulin/pills)
# BPQ: blood pressure questionnaire (self-report + meds)
# SMQ: smoking
# ALQ: alcohol
# PAQ: physical activity
#
# Note: some lab file names differ in older cycles; this script is "best effort":
# it will skip missing files and proceed.


# ---- Helper: safe NHANES download ----
safe_nhanes <- function(table_name) {
  message("Downloading: ", table_name)
  out <- tryCatch(
    nhanesA::nhanes(table_name, translated = FALSE, cleanse_numeric = TRUE),
    error = function(e) {
      message("  -> Skipped (not found or download error): ", table_name)
      NULL
    }
  )
  if (!is.null(out)) out <- dplyr::as_tibble(out)
  out
}


# ---- Download all requested files for all cycles ----
downloaded <- purrr::imap(cycle_suffix, ~{
  suf <- .x
  cyc <- .y
  message("\n=== Cycle ", cyc, " (", suf, ") ===")
  tbls <- tables_for_cycle(suf)
  
  dfs <- purrr::map(tbls, safe_nhanes)
  names(dfs) <- tbls
  
  # Tag cycle
  dfs <- purrr::map(dfs, ~ if (!is.null(.x)) dplyr::mutate(.x, .cycle = suf) else NULL)
  dfs
})

merge_cycle <- function(dfs_named) {
  dfs <- dfs_named[!purrr::map_lgl(dfs_named, is.null)]
  if (length(dfs) == 0) return(NULL)
  
  dfs <- purrr::keep(dfs, ~ "SEQN" %in% names(.x))
  if (length(dfs) == 0) return(NULL)
  
  purrr::reduce(dfs, ~ dplyr::full_join(.x, .y, by = c("SEQN", ".cycle")))
}

merged_by_cycle <- purrr::map(downloaded, merge_cycle)
merged_all <- dplyr::bind_rows(merged_by_cycle)

# ---- Basic cleaning: keep adults (optional) ----
# Many outcomes are adult-focused; change threshold if you want all ages.
nhanes <- merged_all %>%
  # Standardize key demographic variables if present
  mutate(
    age = if ("RIDAGEYR" %in% names(.)) RIDAGEYR else NA_real_,
    sex = if ("RIAGENDR" %in% names(.)) factor(RIAGENDR, levels = c(1,2), labels = c("Male","Female")) else NA,
    # Race/ethnicity coding varies slightly by cycle; RIDRETH1 common for these years
    race_eth = if ("RIDRETH1" %in% names(.)) RIDRETH1 else NA_real_,
    education = if ("DMDEDUC2" %in% names(.)) DMDEDUC2 else NA_real_,
    income_pir = if ("INDFMPIR" %in% names(.)) INDFMPIR else NA_real_
  ) %>%
  # Keep adults by default (you can comment this out)
  filter(is.na(age) | age >= 20)

# ---- Blood pressure: compute mean SBP/DBP from up to 4 readings ----
bp_sys_vars <- c("BPXSY1","BPXSY2","BPXSY3","BPXSY4")
bp_dia_vars <- c("BPXDI1","BPXDI2","BPXDI3","BPXDI4")

existing_sys <- bp_sys_vars[bp_sys_vars %in% names(nhanes)]
existing_dia <- bp_dia_vars[bp_dia_vars %in% names(nhanes)]

nhanes <- nhanes %>%
  rowwise() %>%
  mutate(
    sbp_mean = if (length(existing_sys) > 0) mean(c_across(all_of(existing_sys)), na.rm = TRUE) else NA_real_,
    dbp_mean = if (length(existing_dia) > 0) mean(c_across(all_of(existing_dia)), na.rm = TRUE) else NA_real_
  ) %>%
  ungroup()

harmonize_val <- function(data, ..., default_col = NA_real_) {
  vars <- c(...)
  # 1. Create missing columns as NA so dplyr doesn't error
  for (v in vars) {
    if (!v %in% names(data)) {
      data[[v]] <- NA_real_
    }
  }
  # 2. Coalesce them (take the first non-NA value found)
  # We use standard evaluation with dplyr verbs
  data %>%
    mutate(temp_out = coalesce(!!!syms(vars))) %>%
    pull(temp_out)
}

calc_pa_min <- function(do_var, days_var, min_var) {
  case_when(
    do_var == 2 ~ 0,  # Answered "No" -> 0 minutes
    do_var == 1 & is.na(days_var) ~ NA_real_, # "Yes" but missing days -> NA
    do_var == 1 & is.na(min_var) ~ NA_real_,  # "Yes" but missing mins -> NA
    do_var == 1 ~ days_var * min_var, # Calculation
    TRUE ~ NA_real_ # Missing/Refused "Do" question
  )
}

nhanes <- nhanes %>%
  mutate(
    # --- 1. Sedentary Time (Minutes per day) ---
    # Question: "Minutes sedentary activity" (PAD680)
    # Must clean: 7777 (Refused) and 9999 (Don't Know) -> NA
    sedentary_min = harmonize_val(., "PAD680"),
    sedentary_min = if_else(sedentary_min >= 7777, NA_real_, sedentary_min),
    
    # --- 2. Recreational Activity (Leisure) ---
    # Often the most predictive of health-conscious behavior.
    
    # Vigorous Recreation (e.g., running, basketball)
    # PAQ650 (Do it?), PAQ655 (Days), PAD660 (Minutes)
    rec_vig_min = calc_pa_min(harmonize_val(., "PAQ650"), 
                              harmonize_val(., "PAQ655"), 
                              harmonize_val(., "PAD660")),
    
    # Moderate Recreation (e.g., brisk walk, golf, swimming)
    # PAQ665 (Do it?), PAQ670 (Days), PAD675 (Minutes)
    rec_mod_min = calc_pa_min(harmonize_val(., "PAQ665"), 
                              harmonize_val(., "PAQ670"), 
                              harmonize_val(., "PAD675")),
    
    # Total Recreational MVPA (Minutes/Week)
    pa_rec_total = rec_vig_min + rec_mod_min,
    
    # --- 3. Work & Transport Activity (Optional but thorough) ---
    # Vigorous Work
    work_vig_min = calc_pa_min(harmonize_val(., "PAQ605"), 
                               harmonize_val(., "PAQ610"), 
                               harmonize_val(., "PAD615")),
    # Moderate Work
    work_mod_min = calc_pa_min(harmonize_val(., "PAQ620"), 
                               harmonize_val(., "PAQ625"), 
                               harmonize_val(., "PAD630")),
    # Transportation (Walking/Biking to work/store)
    trans_min = calc_pa_min(harmonize_val(., "PAQ635"), 
                            harmonize_val(., "PAQ640"), 
                            harmonize_val(., "PAD645")),
    
    # Global Physical Activity (Minutes/Week)
    # Summing all domains (Work + Transport + Recreation)
    pa_total_weekly = rec_vig_min + rec_mod_min + work_vig_min + work_mod_min + trans_min,
    # --- Anthropometrics ---
    bmi = harmonize_val(., "BMXBMI"),
    waist_cm = harmonize_val(., "BMXWAIST"),
    
    # --- Demographics ---
    # Poverty Income Ratio (0-5 scale). 
    # <1 is below poverty line. 5 is wealthy (capped).
    pir = harmonize_val(., "INDFMPIR"),
    
    # --- Sleep Duration ---
    # 2005-2014: SLD010H (How much sleep do you get)
    # 2015-2018: SLD012 (Sleep hours)
    sleep_hours = harmonize_val(., "SLD012", "SLD010H"),
    
    # --- Smoking Status (Composite) ---
    # SMQ020: Smoked at least 100 cigarettes in life? (1=Yes, 2=No)
    # SMQ040: Do you now smoke cigarettes? (1=Every day, 2=Some days, 3=Not at all)
    smoke_ever = harmonize_val(., "SMQ020"),
    smoke_now  = harmonize_val(., "SMQ040"),
    
    smoking_status = case_when(
      smoke_now %in% c(1, 2) ~ "Current",
      smoke_ever == 1 & smoke_now == 3 ~ "Former",
      smoke_ever == 2 ~ "Never",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Never", "Former", "Current")),
    
    # --- Liver Enzymes (AST/ALT) ---
    # Usually consistent in BIOPRO as LBXSATSI (ALT) and LBXSASSI (AST)
    # We add LBXALT/LBXAST as backups just in case older/different naming creeps in
    alt = harmonize_val(., "LBXSATSI", "LBXALT"),
    ast = harmonize_val(., "LBXSASSI", "LBXAST"),
    
    # --- CRP (The tricky one) ---
    # 2007-2010: LBXCRP (Standard CRP)
    # 2015-2018: LBXHSCRP (High-Sensitivity CRP)
    # 2011-2014: Data not collected (will return NA)
    crp = harmonize_val(., "LBXHSCRP", "LBXCRP"),
    
    # --- Metabolic / Kidney ---
    # Uric Acid (umol/L) - LBXSUA is usually mg/dL, LBXSUASI is umol/L.
    # Check units carefully. Usually stick to SI (LBXSUASI).
    uric_acid = harmonize_val(., "LBXSUASI", "LBXSUA"),
    
    # Insulin (pmol/L)
    insulin_raw = harmonize_val(., "LBDINSI"),       # Try SI first
    insulin_old = harmonize_val(., "LBXIN"),         # Try conventional
    
    # Final Insulin (Unified to pmol/L)
    insulin = case_when(
      !is.na(insulin_raw) ~ insulin_raw,
      !is.na(insulin_old) ~ insulin_old * 6.00, # Convert old units
      TRUE ~ NA_real_
    ),
    
    # --- Vitamin D (The Complex One) ---
    # NHANES suggests a formula to convert old Vitamin D (2001-2006) to current LC/MS method.
    # For 2007+, usually LBDVIDMS is available directly or requires minor adjustment.
    # We will grab LBDVIDMS (Total 25OH Vit D).
    vitD = harmonize_val(., "LBDVIDMS", "LBXVIDMS", "LBDTCSI"),
    
    # --- Calculated Features for Prediction ---
    # HOMA-IR: (Insulin * Glucose) / 22.5 (if glucose is mmol/L) or / 405 (if mg/dL)
    # NHANES Glucose (LBXGLU) is often mg/dL. 
    # Check units: LBXGLU is mg/dL usually. LBDGLUSI is mmol/L.
    glucose_mmol = harmonize_val(., "LBDGLUSI"),
    glucose_mmol = if_else(is.na(glucose_mmol), harmonize_val(., "LBXGLU") / 18.016, glucose_mmol),
    homa_ir = (insulin * glucose_mmol) / 22.5
  ) 


# ---- Outcome 1: Diabetes ----
# Components (commonly present):
# - HbA1c: LBXGH (GHB file)  :contentReference[oaicite:7]{index=7}
# - Fasting glucose: LBXGLU (GLU file, fasting subsample) :contentReference[oaicite:8]{index=8}
# - Self-reported diabetes: DIQ010 (ever told you have diabetes; standard DIQ variable)
# - Diabetes meds: DIQ050 (taking insulin) and DIQ070 (taking diabetic pills) :contentReference[oaicite:9]{index=9}
#
# Rule: diabetes = 1 if ANY of:
# - HbA1c >= 6.5
# - fasting glucose >= 126 mg/dL
# - self-report diabetes yes
# - taking insulin or diabetic pills yes
#
# Note: fasting glucose is missing for most participants (subsample).
nhanes <- nhanes %>%
  mutate(
    a1c = if ("LBXGH" %in% names(.)) LBXGH else NA_real_,
    fpg = if ("LBXGLU" %in% names(.)) LBXGLU else NA_real_,
    dm_self = if ("DIQ010" %in% names(.)) DIQ010 else NA_real_,
    dm_insulin = if ("DIQ050" %in% names(.)) DIQ050 else NA_real_,
    dm_pills = if ("DIQ070" %in% names(.)) DIQ070 else NA_real_,
    
    diabetes = case_when(
      !is.na(a1c) & a1c >= 6.5 ~ 1L,
      !is.na(fpg) & fpg >= 126 ~ 1L,
      dm_self == 1 ~ 1L,
      dm_insulin == 1 ~ 1L,
      dm_pills == 1 ~ 1L,
      TRUE ~ 0L
    )
  )

# ---- Outcome 2 (alternative): Hypertension ----
# Components:
# - Self-report HTN: BPQ020 (ever told high BP) :contentReference[oaicite:10]{index=10}
# - BP meds: BPQ050A (now taking prescribed medicine) :contentReference[oaicite:11]{index=11}
# - Measured BP: sbp_mean/dbp_mean from BPX exam
#
# Common definition (ACC/AHA 2017 style): SBP>=130 or DBP>=80 or meds or self-report.
# You can tighten to SBP>=140/DBP>=90 if you want a “classic” definition.
nhanes <- nhanes %>%
  mutate(
    htn_self = if ("BPQ020" %in% names(.)) BPQ020 else NA_real_,
    htn_meds = if ("BPQ050A" %in% names(.)) BPQ050A else NA_real_,
    
    hypertension = case_when(
      !is.na(sbp_mean) & sbp_mean >= 130 ~ 1L,
      !is.na(dbp_mean) & dbp_mean >= 80  ~ 1L,
      htn_meds == 1 ~ 1L,
      htn_self == 1 ~ 1L,
      TRUE ~ 0L
    )
  )

mlnhanes <- nhanes %>%
  set_variable_labels(
    # --- Demographics (NHANES Raw) ---
    SEQN      = "Respondent sequence number",
    SDDSRVYR  = "Data release cycle",
    RIDSTATR  = "Interview/Examination status",
    RIAGENDR  = "Gender",
    RIDAGEYR  = "Age in years at screening",
    RIDRETH1  = "Race/Hispanic origin",
    DMDEDUC2  = "Education level (Adults 20+)",
    DMDMARTL  = "Marital status",
    DMDHHSIZ  = "Total number of people in the Household",
    DMDFMSIZ  = "Total number of people in the Family",
    DMDHRGND  = "HH Ref Person Gender",
    SIALANG   = "Language of SP Interview",
    SIAPROXY  = "Proxy used in Interview?",
    SIAINTRP  = "Interpreter used in Interview?",
    
    # --- Weights & Survey Design ---
    WTINT2YR  = "Full sample 2 year interview weight",
    WTMEC2YR  = "Full sample 2 year MEC exam weight",
    SDMVPSU   = "Masked variance pseudo-PSU",
    SDMVSTRA  = "Masked variance pseudo-stratum",
    .cycle    = "NHANES Cycle Year",
    
    # --- Questionnaire Data ---
    DIQ010    = "Doctor told you have diabetes",
    BPQ020    = "Ever told you had high blood pressure",
    SMDUPCA   = "Supplements used (count)", # Assuming this derived count
    SMD100BR  = "Cigarettes brand",
    SMAQUEX2  = "Questionnaire Mode",
    
    # --- Physical Activity (PAQ) ---
    PAQ605    = "Vigorous work activity",
    PAQ620    = "Moderate work activity",
    PAQ635    = "Walk or bicycle for travel",
    PAQ650    = "Vigorous recreational activities",
    PAQ665    = "Moderate recreational activities",
    
    LBXGH = "HbA1c (%)",
    LBXGLU = "Fasting plasma glucose (mg/dL)",
    
    # --- Derived / Renamed Variables ---
    age       = "Age (years)",
    sex       = "Sex",
    race_eth  = "Race and Ethnicity",
    sbp_mean = "Mean systolic blood pressure (mmHg)",
    dbp_mean = "Mean diastolic blood pressure (mmHg)",
    education = "Education Level",
    income_pir= "Ratio of family income to poverty",
    dm_self   = "Self-reported Diabetes",
    diabetes  = "Diabetes Status (Composite / clinical definition))",
    htn_self  = "Self-reported Hypertension",
    hypertension = "Hypertension Status (Composite)")


mlnhanes <- nhanes %>%
  mutate(
    sex = factor(sex, levels = c("Male", "Female")),
    race_eth = factor(race_eth,
                      levels = c(1, 2, 3, 4, 5),
                      labels = c(
                        "Mexican American",
                        "Other Hispanic",
                        "Non-Hispanic White",
                        "Non-Hispanic Black",
                        "Other race"
                      )),
    education = factor(
      education,
      levels = c(1, 2, 3, 4, 5),
      labels = c(
        "<9th grade",
        "9–11th grade",
        "High school",
        "Some college",
        "College graduate"
      )
    ),
    SIALANG=factor(SIALANG,
                   levels=c(1,2),
                   labels=c('English', 'Other')),
  )

mlnhanes <- mlnhanes %>%
  mutate(
    PAQ605 = case_when(
      PAQ605 == 1 ~ "Yes",
      PAQ605 == 2 ~ "No",
      TRUE ~ NA_character_
    ),
    PAQ620 = case_when(
      PAQ620 == 1 ~ "Yes",
      PAQ620 == 2 ~ "No",
      TRUE ~ NA_character_
    ),
    PAQ635 = case_when(
      PAQ635 == 1 ~ "Yes",
      PAQ635 == 2 ~ "No",
      TRUE ~ NA_character_
    ),
    PAQ650 = case_when(
      PAQ650 == 1 ~ "Yes",
      PAQ650 == 2 ~ "No",
      TRUE ~ NA_character_
    ),
    PAQ665 = case_when(
      PAQ665 == 1 ~ "Yes",
      PAQ665 == 2 ~ "No",
      TRUE ~ NA_character_
    ),
    
    PAQ605 = factor(PAQ605, levels = c("No", "Yes")),
    PAQ620 = factor(PAQ620, levels = c("No", "Yes")),
    PAQ635 = factor(PAQ635, levels = c("No", "Yes")),
    PAQ650 = factor(PAQ650, levels = c("No", "Yes")),
    PAQ665 = factor(PAQ665, levels = c("No", "Yes")),
    DMDMARTL = case_when(
      DMDMARTL == 1 ~ "Married",
      DMDMARTL == 2 ~ "Widowed",
      DMDMARTL == 3 ~ "Divorced",
      DMDMARTL == 4 ~ "Separated",
      DMDMARTL == 5 ~ "Never married",
      DMDMARTL == 6 ~ "Living with partner",
      TRUE ~ NA_character_
    ),
    DMDMARTL = factor(
      DMDMARTL,
      levels = c(
        "Married",
        "Living with partner",
        "Never married",
        "Divorced",
        "Separated",
        "Widowed"
      )
    ),
    alcohol = factor(
      case_when(
        ALQ101== 1 ~ "Yes",
        ALQ101 == 2 ~ "No",
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    ),
    SIAPROXY = case_when(
      SIAPROXY == 1 ~ "Yes",
      SIAPROXY == 2 ~ "No",
      TRUE ~ NA_character_
    ),
    SIAPROXY = factor(SIAPROXY, levels = c("No", "Yes")),
    SIAINTRP = case_when(
      SIAINTRP == 1 ~ "Yes",
      SIAINTRP == 2 ~ "No",
      TRUE ~ NA_character_
    ),
    SIAINTRP = factor(SIAINTRP, levels = c("No", "Yes"))
  )

# ---- Optional: quick sanity checks ----
message("\nDiabetes prevalence (unweighted, adults): ", round(mean(nhanes$diabetes, na.rm = TRUE), 3))
message("Hypertension prevalence (unweighted, adults): ", round(mean(nhanes$hypertension, na.rm = TRUE), 3))
message("N rows: ", nrow(nhanes), "  |  N cols: ", ncol(nhanes))

# ---- Save ----
#dir.create("data", showWarnings = FALSE)

message("\nSaved: data/nhanes_2007_2018_ml_ready.rds")

#usethis::use_data(mlnhanes, overwrite = TRUE)
usethis::use_data(mlnhanes, overwrite = TRUE)
saveRDS(mlnhanes, file = file.path("inst", "extdata", "mlnhanes.Rds"))
