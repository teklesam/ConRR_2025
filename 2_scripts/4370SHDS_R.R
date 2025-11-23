
# Samart5 Factory Worker Survival Analysis
# Student ID: XXXXXXX4370
# Date: 2025-11-21
# 
# Research Question: Does occupational exposure to Samart5 increase mortality
# risk in factory workers, and how does this depend on exposure level?


#
# 1. SETUP AND PACKAGE LOADING--------------------------------------------------
#

# Install pacman if not available (for convenient package management)
if (!require("pacman")) install.packages("pacman")

# Load required packages
pacman::p_load(
  rio,          # Import/export data in various formats
  here,         # Manage file paths robustly
  tidyverse,    # Data manipulation and visualization
  forcats,      # Factor manipulation
  skimr,        # Data summaries
  survival,     # Survival analysis
  survminer,    # Enhanced survival plots
  lubridate,    # Date handling
  gtsummary,    # Publication-ready tables
  labelled,     # Variable labels
  Hmisc,        # cut2() for creating quantiles
  patchwork,    # Combine plots
  ggcorrplot,   # Correlation matrices
  car,          # VIF diagnostics
  flextable     # Export tables
)

#
# 2. DATA IMPORT----------------------------------------------------------------
#

# Import dataset
Samart5.data <- import(here('1_data', 'factorydata.csv'))  # This importing technique best works in an R project file.

# Study design details:
# - Cohort study of factory workers
# - Baseline: 31st December 2012
# - Follow-up: 1st January 2013 to 31st December 2015 (3 years)
# - Events: 124 deaths
# - Censored: 98 emigrated, 1825 alive at study end

#
# 3. DATA CLEANING--------------------------------------------------------------
#

#
# 3.1 Calculate follow-up times
#

# Define study period
first.day.ffup <- ymd("2012-12-31")  # Baseline
last.day.ffup <- ymd("2015-12-31")   # End of follow-up

# Add baseline date
Samart5.data$first.day.ffup <- first.day.ffup

# Calculate end of follow-up for each worker
# Workers alive at study end: use last.day.ffup
# Workers who died/emigrated: use baseline + days until event
Samart5.data$last.day.ffup <- if_else(
  is.na(Samart5.data$tim),
  last.day.ffup,
  first.day.ffup + days(Samart5.data$tim)
)

# Calculate follow-up duration in days
Samart5.data$cleaned.time <- as.numeric(
  Samart5.data$last.day.ffup - Samart5.data$first.day.ffup
)

#
# 3.2 Check for missing data
#

# Comprehensive data summary
skim(Samart5.data)

# Result: No missing values in ALL variables

#
# 3.3 Create labeled factors variables - important for legend of plots and tables.
#

Samart5.data <- Samart5.data %>% 
  mutate(
    # Sex labels
    sex.label = case_when(
      sex == '0' ~ 'Male',
      sex == '1' ~ 'Female'
    ),
    
    # Job role labels
    role.label = case_when(
      role == '0' ~ 'Administrative/Managerial',
      role == '1' ~ 'Maintenance',
      role == '2' ~ 'Production'
    ),
    
    # Smoking status labels
    smoke.label = case_when(
      smoke == '0' ~ 'Never smoked',
      smoke == '1' ~ 'Ex-smoker',
      smoke == '2' ~ 'Current smoker'
    ),
    smoke.label = fct_relevel(smoke.label, 'Never smoked', 
                              'Ex-smoker', 'Current smoker'),  # This assumed to make 'Never smoked' category as a reference later, no difference in risk between Ex-smokers and current smokers is assumed.
    
    # Event indicator (0 = censored, 1 = died)
    event = if_else(fail == TRUE, 1, 0),
    
    # Event label for tables
    fail.label = if_else(fail == FALSE, 'Admin censored', 'Died'),
    
    # Composite exposure variable
    # Combines job role intensity with or without above/below median duration of employment
    exposure.cat = case_when(
      role == 0 ~ "Low (Admin)",                                  
      role == 1 & yrswork < median(yrswork[role == 1]) ~                   
        "Medium-Low (Maintenance <10y)",
      role == 1 & yrswork >= median(yrswork[role == 1]) ~ 
        "Medium-High (Maintenance ≥10y)",
      role == 2 & yrswork < median(yrswork[role == 2]) ~ 
        "High (Production <7y)",
      role == 2 & yrswork >= median(yrswork[role == 2]) ~ 
        "Very High (Production ≥7y)"
    ),
    exposure.cat = fct_relevel(exposure.cat,                       # This is turned to ordinal data
                               "Low (Admin)",
                               "Medium-Low (Maintenance <10y)",
                               "Medium-High (Maintenance >=10y)",
                               "High (Production <7y)",
                               "Very High (Production >=7y)"),
    
    # Follow-up time in years
    time.years = cleaned.time / 365.25,                           # This will help to create a year based survival object to generate summaries by years.
    
    # Convert to factors
    sex.label = as.factor(sex.label),
    role.label = as.factor(role.label),
    smoke.label = as.factor(smoke.label),
    fail.label = as.factor(fail.label)
  )

#
# 3.4 Add variable labels for tables
#

Samart5.data <- Samart5.data %>% 
  set_variable_labels(
    sex.label = "Sex",
    age = "Age (in years at 31/12/2012)",
    startage = "Age at Employment (years)",
    yrswork = "Duration of employment (years)",
    role.label = 'Employment role',
    smoke.label = 'Smoking status at 31st December 2012',
    fail.label = 'Died or Censored',
    event = 'Died(1) or Censored (0)',
    cleaned.time = 'Follow-up Time (Days)',
    exposure.cat = 'Samart5 Cumulative Exposure Level',
    time.years = 'Follow-up Time (Years)'
  )

#
# 3.5 Create clean dataset
#

# Remove redundant coded variables, keep labeled versions

Samart5.data.cleaned <- Samart5.data %>%
  select(-c(sex, role, smoke, tim, fail, first.day.ffup, last.day.ffup))

#
# 4. EXPLORATORY DATA ANALYSIS--------------------------------------------------
#

#
# 4.1 Descriptive statistics table
#

# Create table stratified by vital status

summary.table <- Samart5.data.cleaned %>% 
  select(-c(event, cleaned.time)) %>%                                # This table is meant for charaterization of the workers
  tbl_summary(
    by = fail.label,
    statistic = list(
      all_continuous() ~ "mean±sd = {mean} (±{sd}), median[IQR] = {median} ({p25}-{p75})"   # Summarize continous by median, mean and their measures of dispersion
    ),
    digits = all_continuous() ~ 2                                    # round the decimal to 2 digits
  ) %>% 
  add_overall() %>%
  bold_labels() %>%
  add_p() %>%
  bold_p()

summary.table

#
# 4.2 Distribution plots
#
# Examine outlier and distribution shape of the continous variables
# Set up plotting area
par(mfrow = c(1, 3))

# Age distribution
hist(Samart5.data.cleaned$age, 
     breaks = 40, 
     probability = TRUE,
     main = "Age at Baseline",
     xlab = "Age (years)",
     col = "lightblue",
     border = "white")

# Years worked distribution
hist(Samart5.data.cleaned$yrswork, 
     breaks = 40, 
     probability = TRUE,
     main = "Years Worked",
     xlab = "Years",
     col = "lightcoral",
     border = "white")

# Age at employment distribution
hist(Samart5.data.cleaned$startage, 
     breaks = 40, 
     probability = TRUE,
     main = "Age at Employment",
     xlab = "Age (years)",
     col = "lightgreen",
     border = "white")

# Reset plotting area
par(mfrow = c(1, 1))

#
# 4.3 Create tertiles for KM curves
#

# Age tertiles (for checking fitting KM curves and log-log plots for PH assumption testing)
Samart5.data.cleaned <- Samart5.data.cleaned %>%
  mutate(
    age_q = cut2(age, g = 3),
    startage_q = cut2(startage, g = 3),
    yrswork_q = cut2(yrswork, g = 3)
  )

#
# 5. SURVIVAL ANALYSIS - KAPLAN-MEIER CURVES------------------------------------
#

#
# 5.1 Create survival object
#

# Combine time and event indicator
Samart5.survobj <- Surv(
  time = Samart5.data.cleaned$cleaned.time, 
  event = Samart5.data.cleaned$event
)

#
# 5.2 Overall survival curve
#

# Fit Kaplan-Meier for entire cohort
km.overall <- survfit(Samart5.survobj ~ 1, data = Samart5.data.cleaned)

# Create main plot
plot_main <- ggsurvplot(
  km.overall,
  data = Samart5.data.cleaned,
  conf.int = TRUE,
  surv.scale = "percent",
  break.time.by = 200,
  xlab = "Follow-up days",
  ylab = "Survival Probability",
  risk.table = TRUE,
  palette = "#D32F2F",
  conf.int.style = "step",
  ggtheme = theme_light(base_size = 20) +
    theme(
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 18),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 20)
    )
)$plot

# Create zoomed inset (90-100% survival)
plot_zoom <- ggsurvplot(
  km.overall,
  data = Samart5.data.cleaned,
  conf.int = TRUE,
  surv.scale = "percent",
  break.time.by = 200,
  xlab = "",
  ylab = "",
  risk.table = FALSE,
  palette = "#D32F2F",
  conf.int.style = "step",
  ggtheme = theme_light(base_size = 8)
)$plot +
  coord_cartesian(ylim = c(0.90, 1.00)) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_blank()
  )

# Combine main plot with inset
plot_final <- plot_main +
  inset_element(
    plot_zoom,
    left = 0.45, bottom = 0.18,
    right = 0.98, top = 0.70
  )

# Save figure
ggsave(
  "KM_curve_overall.pdf",
  plot_final,
  width = 10, height = 7, dpi = 300
)

#
# 5.3 Survival by Samart 5 exposure level
#

# Fit Kaplan-Meier stratified by exposure
km.exposure <- survfit(Samart5.survobj ~ exposure.cat, 
                       data = Samart5.data.cleaned)

# Create survival plot
plot.exposure <- ggsurvplot(
  km.exposure,
  data = Samart5.data.cleaned,
  risk.table = FALSE,
  pval = TRUE,
  conf.int = TRUE,
  surv.scale = "percent",
  palette = c("#2E7D32", "#66BB6A", "#2E9FDF", "#E7B800", "#D32F2F"),
  xlab = "Time in days",
  break.time.by = 100,
  conf.int.style = "step",
  surv.median.line = "hv",
  legend.labs = c("Low (Admin)", 
                  "Medium-Low (Maintenance <10y)", 
                  "Medium-High (Maintenance >=10y)", 
                  "High (Production <7y)", 
                  "Very High (Production ≥7y)"),
  ggtheme = theme_light(base_size = 20) +
    theme(
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 18),
      legend.text = element_text(size = 22),
      legend.title = element_text(size = 20)
    )
)

# Zoom to show differences clearly                   # the puporse of zooming the plots isnot to accentuate the difference, but to see the difference in survival curves better.
plot.exposure$plot <- plot.exposure$plot +
  coord_cartesian(ylim = c(0.70, 1.00))

# Save figure
pdf(
  file = here("3_output", "km_exposure_2.pdf"),
  width = 21, height = 15
)
print(plot.exposure)
dev.off()

# Formal log-rank test
survdiff(Samart5.survobj ~ exposure.cat, data = Samart5.data.cleaned)

#
# 5.4 Survival probability table
#

# Estimate survival at 1st, 2nd, and 3rd years                 # As we cant measure the median survival times due to high survial probability in the cohort, we will compare the survival probabilty at 1st, 2nd and 3rd year
survival_table <- Samart5.data.cleaned %>%
  select(time.years, event, exposure.cat, role.label, 
         sex.label, smoke.label, yrswork_q, age_q) %>%
  tbl_survfit(
    y = Surv(time.years, event),
    times = c(1, 2, 3),
    label_header = "**{time} Year Survival (95% CI)**"
  ) %>%
  add_p(test = "logrank") %>%
  add_n() %>%
  bold_p() %>%
  bold_labels()

survival_table

#
# 6. COX PROPORTIONAL HAZARDS MODELS--------------------------------------------
#

#
# 6.1 Univariable models
#

# Exposure
cox_uni_exposure <- coxph(Samart5.survobj ~ exposure.cat, 
                          data = Samart5.data.cleaned)

# Age
cox_uni_age <- coxph(Samart5.survobj ~ age, 
                     data = Samart5.data.cleaned)

# Sex
cox_uni_sex <- coxph(Samart5.survobj ~ sex.label, 
                     data = Samart5.data.cleaned)

# Smoking
cox_uni_smoke <- coxph(Samart5.survobj ~ smoke.label, 
                       data = Samart5.data.cleaned)

#
# 6.2 Multivariable model
#

# Fit adjusted Cox model
# Note: Starting age excluded due to collinearity (age = startage + yrswork)
cox_multi <- coxph(
  Samart5.survobj ~ exposure.cat + age + sex.label + smoke.label,
  data = Samart5.data.cleaned
)

# Model summary
summary(cox_multi)

#
# 6.3 Check proportional hazards assumption
#

# Schoenfeld residuals test
ph_test <- cox.zph(cox_multi)
print(ph_test)

# Visual inspection
ggcoxzph(ph_test)

# Complementary log-log plots for Samart5 exposure, age, sex, smoking, 
fit_exposure <- survfit(Surv(cleaned.time, event) ~ exposure.cat, 
                        data = Samart5.data.cleaned)

ggsurvplot(
  fit_exposure,
  data = Samart5.data.cleaned,
  fun = "cloglog",
  xlab = "Log(Follow-up time in days)",
  ylab = "log(-log(S(t)))",
  conf.int = FALSE,
  palette = "Dark2",
  ggtheme = theme_light(),
  risk.table = FALSE,
  title = "Proportional Hazards Check: Exposure Level"
)

fit.age_q <- survfit(Surv(cleaned.time, event) ~ age_q, 
                     data = Samart5.data.cleaned)

ggsurvplot(
  fit.age_q,
  data = Samart5.data.cleaned,
  fun = "cloglog",
  xlab = "Log(Follow-up time in days)",
  ylab = "log(-log(S(t)))",
  conf.int = FALSE,
  palette = "Dark2",
  ggtheme = theme_light(),
  legend.title = "Age Tertiles",
  risk.table = FALSE,
  title = "Proportional Hazards Check: Age"
)

fit.sex <- survfit(Surv(cleaned.time, event) ~ sex.label, 
                   data = Samart5.data.cleaned)

ggsurvplot(
  fit.sex,
  data = Samart5.data.cleaned,
  fun = "cloglog",
  xlab = "Log(Follow-up time in days)",
  ylab = "log(-log(S(t)))",
  conf.int = FALSE,
  palette = "Dark2",
  ggtheme = theme_light(),
  legend.title = "Sex",
  risk.table = FALSE,
  title = "Proportional Hazards Check: Sex"
)

fit.startage_q <- survfit(Surv(cleaned.time, event) ~ startage_q, 
                          data = Samart5.data.cleaned)

ggsurvplot(
  fit.startage_q,
  data = Samart5.data.cleaned,
  fun = "cloglog",
  xlab = "Log(Follow-up time in days)",
  ylab = "log(-log(S(t)))",
  conf.int = FALSE,
  palette = "Dark2",
  ggtheme = theme_light(),
  legend.title = "Starting Age Tertiles",
  risk.table = FALSE,
  title = "Proportional Hazards Check: Starting Age"
)

fit.role <- survfit(Surv(cleaned.time, event) ~ role.label, 
                    data = Samart5.data.cleaned)

ggsurvplot(
  fit.role,
  data = Samart5.data.cleaned,
  fun = "cloglog",
  xlab = "Log(Follow-up time in days)",
  ylab = "log(-log(S(t)))",
  conf.int = FALSE,
  palette = "Dark2",
  ggtheme = theme_light(),
  legend.title = "Employment Role",
  risk.table = FALSE,
  title = "Proportional Hazards Check: Employment Role"
)

fit.smoke <- survfit(Surv(cleaned.time, event) ~ smoke.label, 
                     data = Samart5.data.cleaned)

ggsurvplot(
  fit.smoke,
  data = Samart5.data.cleaned,
  fun = "cloglog",
  xlab = "Log(Follow-up time in days)",
  ylab = "log(-log(S(t)))",
  conf.int = FALSE,
  palette = "Dark2",
  ggtheme = theme_light(),
  legend.title = "Smoking Status",
  risk.table = FALSE,
  title = "Proportional Hazards Check: Smoking Status"
)

#
# 6.4 Model diagnostics
#

# Check multicollinearity
vif(cox_multi)

#
# 7. CREATE RESULTS TABLES------------------------------------------------------
#

#
# 7.1 Univariable results table
#

tbl_uni <- Samart5.data.cleaned %>%
  select(exposure.cat, age, sex.label, smoke.label) %>%
  tbl_uvregression(                                        # Dedicated for univariate table creation
    method = coxph,
    y = Samart5.survobj,
    exponentiate = TRUE,                                   # Present the coffecients after exponentiatning them
    label = list(
      exposure.cat ~ "Samart5 Cumulative Exposure Level",
      age ~ "Age (in years at 31/12/2012)",
      sex.label ~ "Sex",
      smoke.label ~ "Smoking status at 31st December 2012"
    )
  ) %>%
  add_global_p() %>%
  bold_p() %>%
  bold_labels()

#
# 7.2 Multivariable results table
#

tbl_multi <- tbl_regression(
  cox_multi,
  exponentiate = TRUE,
  label = list(
    exposure.cat ~ "Samart5 Cumulative Exposure Level",
    age ~ "Age (in years at 31/12/2012)",
    sex.label ~ "Sex",
    smoke.label ~ "Smoking status at 31st December 2012"
  )
) %>%
  add_global_p() %>%                           #NB. The global p value is added, because it  tests the null hypothesis that all coefficients for that multi-level categorical variable are equal to zero. 95% CI also does hypothesis testing, its provided in the table for every cofficient
  add_glance_table(include = AIC) %>%
  add_vif() %>%
  bold_p() %>%
  bold_labels()

#
# 7.3 Combined table
#

final_table <- tbl_merge(
  tbls = list(tbl_uni, tbl_multi),
  tab_spanner = c("**Unadjusted**", "**Adjusted**")
) %>%
  modify_caption("**Univariate and Multivariable Cox Models**")

final_table

# Export to Word  to be edited in overleaf
# final_table %>%
#   as_flex_table() %>%
#   save_as_docx(path = here("3_output", "cox_models_table.docx"))

#
# 8. SESSION INFORMATION-------------
#

# Document R version and packages for reproducibility
sessionInfo()
> sessionInfo()
# R version 4.5.1 (2025-06-13 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26100)
# 
# Matrix products: default
# LAPACK version 3.12.1
# 
# locale:
#   [1] LC_COLLATE=English_United Kingdom.utf8  LC_CTYPE=C                              LC_MONETARY=English_United Kingdom.utf8
# [4] LC_NUMERIC=C                            LC_TIME=English_United Kingdom.utf8    
# system code page: 65001
# 
# time zone: Europe/London
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] flextable_0.9.10   car_3.1-3          carData_3.0-5      ggcorrplot_0.1.4.1 patchwork_1.3.2    Hmisc_5.2-4       
# [7] labelled_2.16.0    gtsummary_2.4.0    survminer_0.5.1    ggpubr_0.6.2       survival_3.8-3     skimr_2.2.1       
# [13] lubridate_1.9.4    forcats_1.0.1      stringr_1.5.2      dplyr_1.1.4        purrr_1.1.0        readr_2.1.5       
# [19] tidyr_1.3.1        tibble_3.3.0       ggplot2_4.0.0      tidyverse_2.0.0    here_1.0.2         rio_1.2.4         
# [25] pacman_0.5.1      
# 
# loaded via a namespace (and not attached):
#   [1] gridExtra_2.3           rlang_1.1.6             magrittr_2.0.4          compiler_4.5.1          systemfonts_1.3.1      
# [6] vctrs_0.6.5             reshape2_1.4.5          pkgconfig_2.0.3         crayon_1.5.3            fastmap_1.2.0          
# [11] backports_1.5.0         labeling_0.4.3          utf8_1.2.6              KMsurv_0.1-6            rmarkdown_2.30         
# [16] markdown_2.0            tzdb_0.5.0              haven_2.5.5             ragg_1.5.0              xfun_0.54              
# [21] litedown_0.8            jsonlite_2.0.0          uuid_1.2-1              broom_1.0.10            cluster_2.1.8.1        
# [26] R6_2.6.1                stringi_1.8.7           RColorBrewer_1.1-3      rpart_4.1.24            Rcpp_1.1.0             
# [31] knitr_1.50              zoo_1.8-14              base64enc_0.1-3         parameters_0.28.2       R.utils_2.13.0         
# [36] Matrix_1.7-3            splines_4.5.1           nnet_7.3-20             timechange_0.3.0        tidyselect_1.2.1       
# [41] rstudioapi_0.17.1       abind_1.4-8             yaml_2.3.10             ggtext_0.1.2            lattice_0.22-7         
# [46] plyr_1.8.9              withr_3.0.2             bayestestR_0.17.0       S7_0.2.0                askpass_1.2.1          
# [51] evaluate_1.0.5          foreign_0.8-90          zip_2.3.3               xml2_1.4.0              survMisc_0.5.6         
# [56] pillar_1.11.1           rsconnect_1.5.1         checkmate_2.3.3         insight_1.4.2           generics_0.1.4         
# [61] rprojroot_2.1.1         hms_1.1.3               commonmark_2.0.0        scales_1.4.0            xtable_1.8-4           
# [66] glue_1.8.0              gdtools_0.4.4           tools_4.5.1             data.table_1.17.8       ggsignif_0.6.4         
# [71] fs_1.6.6                cowplot_1.2.0           grid_4.5.1              datawizard_1.3.0        cards_0.7.0            
# [76] colorspace_2.1-2        repr_1.1.7              cardx_0.3.0             htmlTable_2.4.3         Formula_1.2-5          
# [81] cli_3.6.5               km.ci_0.5-6             textshaping_1.0.4       officer_0.7.0           fontBitstreamVera_0.1.1
# [86] broom.helpers_1.22.0    gt_1.1.0                gtable_0.3.6            R.methodsS3_1.8.2       rstatix_0.7.3          
# [91] fontquiver_0.2.1        sass_0.4.10             digest_0.6.37           htmlwidgets_1.6.4       farver_2.1.2           
# [96] htmltools_0.5.8.1       R.oo_1.27.1             lifecycle_1.0.4         openssl_2.3.4           fontLiberation_0.1.0   
# [101] gridtext_0.1.5

# END OF ANALYSIS
