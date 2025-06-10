#------------------------------
# MAIN ANALYSIS
#------------------------------

# This script contains code for the main analysis.
# Two different types of models are fitted to the pre-processed data.

#-----------HEADER-------------
# load packages
library(knitr)        # output table editing etc.
library(tidyverse)    # data management and ggplot2
library(mgcv)         # Additive regression models
library(Hmisc)
library(emmeans)
library(dplyr)
library(sjPlot)
theme_set(theme_bw()) # set global theme for ggplot2 plots
source("plot_predictions.R")

# define paths
path_prepData  <- "0_data/2_preparedData/"
path_resTables <- "1_resultFiles/tables/"
path_resImages <- "1_resultFiles/images/"

dat <- read.csv2(paste0(path_prepData, "data_prepared.csv"), fileEncoding = "utf-8")
#------------------------------

#------------------------------ DATA PREPARATION
# ensure that Participant and Topic are encoded as factor variables
dat <- dat %>%
  mutate(Participant = factor(Participant),
         Topic       = factor(Topic))

# MODEL 1 SPECIFICATION

model_details <- T
# response variables
response_vec <- c("TRT","FRT","RRT")
# define for which models the 'letters per word' effect should be estimated
# non-linearly (default: all)
nonlinear_effect <- c("TRT","FRT","RRT") # set to c() to estimate all 'letters per word' effects as linear effects
no_effect <- c() # set to c("TRT") to remove the 'letters per word' effect in TRT-model
# AOI.conditions per model
AOIcond_list <- list(
  "Model1" = c("E_LM_1", "E_LM_2", "E_LM_3", "E_LM_4", "E_LM_5", "E_LM_6",
               "E_PM_1", "E_PM_2", "E_PM_3", "E_PM_4", "E_PM_5", "E_PM_6",
               "H_LM_1", "H_LM_2", "H_LM_3", "H_LM_4", "H_LM_5", "H_LM_6",
               "H_PM_1", "H_PM_2", "H_PM_3", "H_PM_4", "H_PM_5", "H_PM_6",
               "N_LM_1", "N_LM_2", "N_LM_3", "N_LM_4", "N_LM_5", "N_LM_6",
               "N_PM_1", "N_PM_2", "N_PM_3", "N_PM_4", "N_PM_5", "N_PM_6")
)

AOIcond_list_table_3 <- list(
  "Model1" = c("N_LM_1", "N_PM_1", "E_LM_1", "E_PM_1", "H_LM_1", "H_PM_1",
               "N_LM_2", "N_PM_2", "E_LM_2", "E_PM_2", "H_LM_2","H_PM_2",
               "N_LM_3", "N_PM_3", "E_LM_3", "E_PM_3", "H_LM_3", "H_PM_3")
)

labels_customized <- c("N_LM_1" = "L1, lex, full", "N_PM_1" = "L1, pron, full",
                      "E_LM_1" = "L2, lex, full", "E_PM_1" = "L2, pron, full",
                      "H_LM_1" = "LH, lex, full", "H_PM_1" = "LH, pron, full",
                      "N_LM_2" = "L1, lex, relat", "N_PM_2" = "L1, pron, relat",
                      "E_LM_2" = "L2, lex, relat", "E_PM_2" = "L2, pron, relat",
                      "H_LM_2" = "LH, lex, relat", "H_PM_2" = "LH, pron, relat",
                      "N_LM_3" = "L1, lex, refer", "N_PM_3" = "L1, pron, refer",
                      "E_LM_3" = "L2, lex, refer", "E_PM_3" = "L2, pron, refer",
                      "H_LM_3" = "LH, lex, refer", "H_PM_3" = "LH, pron, refer")


AOIcond_vec <- sort(unique(unlist(AOIcond_list)))
dat <- dat %>% filter(AOI.condition %in% AOIcond_vec) %>%
  mutate(AOI.condition = factor(AOI.condition, levels = AOIcond_vec))

# calculate mean nLetters.WD
nLetters.WD <- round(mean(dat$nLetters.WD), 2)

#------------------------------ DESCRIPTIVES
# Observations per AOI.condition
freq <- sapply(AOIcond_vec, function(AOIcond) {
  length(unique(dat$observation[dat$AOI.condition == AOIcond]))
})
tab <- data.frame("AOI.condition" = AOIcond_vec,
                  "Freq" = unname(freq))
# kable(tab)
#------------------------------

#-----------MODEL 1------------- 

#------------------------------ MODEL 1 ESTIMATION
model_list <- lapply(seq_len(length(AOIcond_list)), function(i) {
  m_list <- lapply(seq_len(length(response_vec)), function(j) {
    message(paste0("Estimate model ",i,"_",response_vec[j],"..."))
    dat_i <- dat %>% mutate(AOI.condition = relevel(AOI.condition, ref = AOIcond_list[[i]][1]))
    fm_covars <- "~ AOI.condition"
    if (!(response_vec[j] %in% no_effect)) { # add 'letters per word' effect
      if (response_vec[j] %in% nonlinear_effect) {
        fm_covars <- paste(fm_covars, "+ s(nLetters.WD, bs='ps', k = 4)")
      } else { 
        fm_covars <- paste(fm_covars, "+ nLetters.WD")
      }
    }
    fm_covars <- paste(fm_covars, "+ s(Topic, bs='re') + s(Participant, bs='re')")
    fm <- as.formula(paste0(response_vec[j], ".WD", fm_covars))
    bam(fm, data = dat_i, method = "REML", nthreads = 10)
  })
  names(m_list) <- response_vec
  m_list
})

names(model_list) <- paste0("Model", seq_len(length(model_list)))
#------------------------------

#------------------------------ MODEL 1 RESULTS
for (i in seq_len(length(AOIcond_list))) {
  # retrieve the y limits for plotting over all response variables
  max_fitted <- sapply(model_list[[i]], function(model) { max(model$fitted.values[model$model$AOI.condition %in% AOIcond_list[[i]]]) })
  ylim <- c(0, max(max_fitted))
  
  for (j in seq_len(length(response_vec))) {
    res_list <- plot_predictions(model = model_list[[i]][[j]],
                                 nLetters.WD = nLetters.WD,
                                 groups_vec = AOIcond_list_table_3[[i]],
                                 model_name = paste("Model",i),
                                 ylim = ylim,
                                 xlab_size = 9,
                                 labels_custom_yes = TRUE,
                                 labels_custom = labels_customized)
    
    # 1) Plot of model predictions
    ggsave(plot = res_list$plot, filename = paste0(path_resImages, "plot_m",i,"_", response_vec[j], ".pdf"), width = 8, height = 5)
    
    # 2) Table of estimates and model predictions
    write.csv(res_list$table, file = paste0(path_resTables, "tab_m",i,"_",response_vec[j],".csv"), row.names = TRUE)
  }
}

# results of model 1 for inspection:
res_mod1_TRT <- read.csv(paste0(path_resTables, "tab_m1_TRT.csv"))
res_mod1_FRT <- read.csv(paste0(path_resTables, "tab_m1_FRT.csv"))
res_mod1_RRT <- read.csv(paste0(path_resTables, "tab_m1_RRT.csv"))
#------------------------------

#------------------------------ VISUALIZATIONS

# plot predictions
dat_pred_plot <- res_mod1_TRT %>%
  mutate(
    speaker = str_extract(X, "^[A-Z]"),
    condition = str_extract(X, "LM|PM"),
    aoi = as.numeric(str_extract(X, "\\d+$"))
  )

vec_titles <- c("Full sentence", "Relational expression", "Referential expression")
plot_by_item <- function(AOI) {
  p <- ggplot(dat_pred_plot %>% filter(aoi == AOI), aes(x = condition, y = TRT.Pred, color = speaker)) +
    geom_point(position = position_dodge(width = 0.15), size = 3) +
    geom_errorbar(aes(ymin = TRT.ci_lower, ymax = TRT.ci_upper),
                  width = 0.2, position = position_dodge(width = 0.15)) +
    geom_line(aes(group = speaker), position = position_dodge(width = 0.15)) +
    scale_color_manual(labels = c("N" = "L1", "E" = "L2", "H" = "LH"),
                       values = c("N" = "#440154", "E" =  "orange", "H" = "#21908C")) +
    labs(title = vec_titles[AOI], y = "predicted TRT per word\nfor fixed number of letters", x = "meaning") +
    scale_x_discrete(labels = c("LM" = "lex", "PM" = "pron")) +
    theme_minimal()
  
  ggsave(plot = p, filename = paste0(path_resImages, "plot_mod1_TRT_AOI", AOI, ".pdf"), 
         width = 5, height = 3)
}

plot_by_item(AOI = 1)
plot_by_item(AOI = 2)
plot_by_item(AOI = 3)

#------------------------------

#-----------MODEL 2------------ 

# DATA PREPARATION
# ----------------
dat$speaker <- factor(substr(dat$condition, 1, 1))
dat$meaning <- factor(substr(dat$AOI.condition, 3, 4))
dat$AOI <- factor(substr(dat$AOI.condition, 6, 6))

# factors with desired level order
dat$speaker <- relevel(dat$speaker, ref = "N")

# MODEL ESTIMATION
# ----------------
mod2_fit <- function(AOI, resp) {
  dat_mod <- dat[dat$AOI == AOI,]
  dat_mod[, paste0(resp, ".L")] <- dat_mod[, resp] / dat_mod$nLetters
  # model estimation
  message(paste0("Estimate model 2_", resp, "_", AOI, "..."))
  form <- as.formula(paste0(resp, ".L ~ speaker * meaning + ",
                            "s(Participant, bs = 're') + s(Topic, bs = 're')"))
  mod <- gam(data = dat_mod, formula = form)
  
  return(mod)
}

mod2_res <- function(AOI, resp) {
  mod <- mod2_fit(AOI, resp)
  assign(paste0("mod2_", resp, "_", AOI), mod, envir = .GlobalEnv)
  
  # res <- round(anova.gam(mod)$pTerms.table, 4)
  # assign(paste0("res_mod2_", resp, "_", AOI), res, envir = .GlobalEnv)
  
  # marginal effects
  contrast_matrix <- list(
    "N_LM vs N_PM" = c(1, 0, 0, -1, 0, 0),
    "E_LM vs E_PM" = c(0, 1, 0, 0, -1, 0),
    "H_LM vs H_PM" = c(0, 0, 1, 0, 0, -1),
    "N_LM vs E_LM" = c(1, -1, 0, 0, 0, 0),
    "N_LM vs H_LM" = c(1, 0, -1, 0, 0, 0),
    "E_LM vs H_LM" = c(0, 1, -1, 0, 0, 0),
    "N_PM vs E_PM" = c(0, 0, 0, 1, -1, 0),
    "N_PM vs H_PM" = c(0, 0, 0, 1, 0, -1),
    "E_PM vs H_PM" = c(0, 0, 0, 0, 1, -1)
  )
  
  emm_contrast <- contrast(emmeans(mod, specs = c("speaker", "meaning")),
                           method = contrast_matrix)
  
  assign(paste0("contrast_mod2_", resp, "_", AOI),
         as.data.frame(emm_contrast), envir = .GlobalEnv)
}

# AOI 1
mod2_res(1, "TRT")
summary(mod2_TRT_1)
kable(contrast_mod2_TRT_1)

mod2_res(1, "FRT")
summary(mod2_FRT_1)
kable(contrast_mod2_FRT_1)

mod2_res(1, "RRT")
summary(mod2_RRT_1)
kable(contrast_mod2_RRT_1)

# AOI 2
mod2_res(2, "TRT")
summary(mod2_TRT_2)
kable(contrast_mod2_TRT_2)

mod2_res(2, "FRT")
summary(mod2_FRT_2)
kable(contrast_mod2_FRT_2)

mod2_res(2, "RRT")
summary(mod2_RRT_2)
kable(contrast_mod2_RRT_2)

# AOI 3
mod2_res(3, "TRT")
summary(mod2_TRT_3)
kable(contrast_mod2_TRT_3)

mod2_res(3, "FRT")
summary(mod2_FRT_3)
kable(contrast_mod2_FRT_3)

mod2_res(3, "RRT")
summary(mod2_RRT_3)
kable(contrast_mod2_RRT_3)

#------------------------------

#------------------------------ VISUALIZATIONS

mod2_effect_plot <- function(AOI, resp) {
  dat_mod <- dat[dat$AOI == AOI,]
  dat_mod[, paste0(resp, ".L")] <- dat_mod[, resp] / dat_mod$nLetters
  
  p <- dat_mod %>% 
    ggplot() +
    aes(x = meaning, color = speaker, group = speaker, y = get(paste0(resp, ".L"))) +
    stat_summary(fun = mean, geom = "point") +
    stat_summary(fun = mean, geom = "line") +
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.1) +
    ylab(paste0("estimates for ", resp, " per letter")) +
    scale_color_manual(labels = c("N" = "L1", "E" = "L2", "H" = "LH"),
                       values = c("N" = "#440154", "E" =  "orange", "H" = "#21908C")) +
    scale_x_discrete(labels = c("LM" = "lex", "PM" = "pron"))
  
  ggsave(plot = p, filename = paste0(path_resImages, "plot_mod2_", resp, "_AOI", AOI, ".pdf"), 
         width = 5, height = 3)
}

mod2_effect_plot(1, "TRT")
mod2_effect_plot(2, "TRT")
mod2_effect_plot(3, "TRT")

#------------------------------
