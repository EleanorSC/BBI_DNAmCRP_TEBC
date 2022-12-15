## ---------------------------

setwd("/Users/eleanorc_worklaptop/repos/BBI_DNAmCRP")

## ----------------------------# 

## load datasets --------------# 
DATA <- read.csv("TEBC_DNAm_data_variables.csv")
PC1_DATA <- read.csv("Eleanor_desc.csv")

## ----------------------------# 

## ----------EXAMINE CORRELATION BETWEEN PC1 gestational age and DNAm CRP-----------# 

PC1 <- PC1_DATA %>% select(ID, PC1)
DATA_FULL <- merge(DATA, PC1,
                   by = "ID")
res <- cor.test(DATA_FULL$scores, DATA_FULL$PC1, 
                method = "pearson")
corr_df <- DATA_FULL %>% select(scores, gest_age)
M <-cor(corr_df, use = "pairwise.complete.obs")
M
## ----------EXAMINE CORRELATION BETWEEN PC1 gestational age and DNAm CRP-----------# 

skimr::skim(DATA_FULL$sex)

DF_214 <- DATA_FULL

DF_214  <-
  DF_214  %>% mutate(
    infant_sex =
      case_when(
        sex == 2 ~ 1,
        sex == 1 ~ 0,
        TRUE ~ 3))


brain_volumes_preterms <- DF_214 %>% filter(Preterm == 1)
brain_volumes_terms <- DF_214 %>% filter(Preterm == 0)

### create a loop

neuroimaging_list <- list(
  "Cortical_grey_matter",
  "Deep_grey_matter",
  "White_matter",
  "Hippocampi_and_Amygdala",
  "Cerebellum",
  "CSF",
  "Ventricles",
  "Brainstem",
  "PSFA",
  "PSMD"
)

neuroimaging_regressions <- list()

### H1 

### Model 1: standard model
i <- 1
for (brain_var in neuroimaging_list) { 
  linear_mod_output <- summary(lm(scale(brain_volumes_preterms[brain_var]) ~ 
                                    scale(gest_age) 
                                  + infant_sex 
                                  + scale(BW_Z) 
                                  + scale(gest_scan) 
                                  + Phase
                                  + scale(scores), 
                                  data = brain_volumes_preterms))
  rownames(linear_mod_output$coefficients)[7] <- brain_var
  neuroimaging_regressions[[i]] <- linear_mod_output
  i <- (i + 1)
}

########

#MODEL 1 (term : age  + infant_sex) This function outputs my various regression summaries into a table
summary_table_h1_preterm <- function(neuroimaging_regressions, neuroimaging_list, brain_var) {
  
  beta <- NULL
  se <- NULL
  pvals <- NULL
  r2 <- NULL
  model <- NULL
  exposure <- NULL
  natal <- NULL
  Preterm <- NULL
  
  for (regression in neuroimaging_regressions) {
    beta <- c(beta, regression$coefficients[7,1])
    se <- c(se, regression$coefficients[7,2])
    pvals <- c(pvals, regression$coefficients[7,4])
    r2 <- c(r2, regression$r.squared)
    #exposure <- c()
  }
  
  return(data.frame(beta = beta,
                    se = se,
                    pvals = pvals,
                    r2 = r2,
                    model = "gest_age  + infant_sex + BW_Z + gest_scan",
                    exposure = "fully_adj",
                    natal = "postnatal",
                    Preterm = "both",
                    brain_var = unlist(neuroimaging_list)
                    
  ))
}


#Here tell the function which regressions to use
H1 <- summary_table_h1_preterm(neuroimaging_regressions, neuroimaging_list, neuroimaging_list)
H1


### H3 
### Model 3: Maternal smoking in pregnancy
i <- 1

for (brain_var in neuroimaging_list) { 
  linear_mod_output <- summary(lm(scale(brain_volumes_preterms[brain_var]) ~ 
                                    scale(gest_age) 
                                  + infant_sex 
                                  + scale(BW_Z) 
                                  + scale(gest_scan) 
                                  + Phase
                                  + Smoked_in_pregnancy
                                  + scale(scores), 
                                  data = brain_volumes_preterms))
  rownames(linear_mod_output$coefficients)[8] <- brain_var
  neuroimaging_regressions[[i]] <- linear_mod_output
  i <- (i + 1)
}

########

#MODEL 1 (term : age  + infant_sex) This function outputs my various regression summaries into a table
summary_table_h2_preterm <- function(neuroimaging_regressions, neuroimaging_list, brain_var) {
  
  beta <- NULL
  se <- NULL
  pvals <- NULL
  r2 <- NULL
  model <- NULL
  exposure <- NULL
  natal <- NULL
  Preterm <- NULL
  
  for (regression in neuroimaging_regressions) {
    beta <- c(beta, regression$coefficients[8,1])
    se <- c(se, regression$coefficients[8,2])
    pvals <- c(pvals, regression$coefficients[8,4])
    r2 <- c(r2, regression$r.squared)
    #exposure <- c()
  }
  
  return(data.frame(beta = beta,
                    se = se,
                    pvals = pvals,
                    r2 = r2,
                    model = "Smoked during pregnancy",
                    exposure = "antenatal",
                    natal = "antenatal",
                    Preterm = "Preterm",
                    brain_var = unlist(neuroimaging_list)
                    
  ))
}


#Here tell the function which regressions to use
H3_smoked_pregnancy <- summary_table_h2_preterm(neuroimaging_regressions, neuroimaging_list, neuroimaging_list)
H3_smoked_pregnancy



### H8 
### Model 8: preeclampsia
i <- 1

for (brain_var in neuroimaging_list) { 
  linear_mod_output <- summary(lm(scale(brain_volumes_preterms[brain_var]) 
                                  ~ scale(gest_age) 
                                  + infant_sex 
                                  + scale(BW_Z) 
                                  + scale(gest_scan) 
                                  + Phase
                                  + preeclampsia
                                  + scale(scores), 
                                  data = brain_volumes_preterms))
  rownames(linear_mod_output$coefficients)[8] <- brain_var
  neuroimaging_regressions[[i]] <- linear_mod_output
  i <- (i + 1)
}

########

#MODEL 1 (term : age  + infant_sex) This function outputs my various regression summaries into a table
summary_table_h2_preterm <- function(neuroimaging_regressions, neuroimaging_list, brain_var) {
  
  beta <- NULL
  se <- NULL
  pvals <- NULL
  r2 <- NULL
  model <- NULL
  exposure <- NULL
  natal <- NULL
  Preterm <- NULL
  
  for (regression in neuroimaging_regressions) {
    beta <- c(beta, regression$coefficients[8,1])
    se <- c(se, regression$coefficients[8,2])
    pvals <- c(pvals, regression$coefficients[8,4])
    r2 <- c(r2, regression$r.squared)
    #exposure <- c()
  }
  
  return(data.frame(beta = beta,
                    se = se,
                    pvals = pvals,
                    r2 = r2,
                    model = "preeclampsia",
                    exposure = "perinatal",
                    natal = "perinatal",
                    Preterm = "Preterm",
                    brain_var = unlist(neuroimaging_list)
                    
  ))
}


#Here tell the function which regressions to use
H8_preeclampsia <- summary_table_h2_preterm(neuroimaging_regressions, neuroimaging_list, neuroimaging_list)
H8_preeclampsia

### H9
### Model 9: HCA
i <- 1

for (brain_var in neuroimaging_list) { 
  linear_mod_output <- summary(lm(scale(brain_volumes_preterms[brain_var]) ~ scale(gest_age) 
                                  + infant_sex 
                                  + scale(BW_Z) 
                                  + scale(gest_scan) 
                                  + Phase
                                  + HCA
                                  + scale(scores), 
                                  data = brain_volumes_preterms))
  rownames(linear_mod_output$coefficients)[8] <- brain_var
  neuroimaging_regressions[[i]] <- linear_mod_output
  i <- (i + 1)
}

########

#MODEL 1 (term : age  + infant_sex) This function outputs my various regression summaries into a table
summary_table_h2_preterm <- function(neuroimaging_regressions, neuroimaging_list, brain_var) {
  
  beta <- NULL
  se <- NULL
  pvals <- NULL
  r2 <- NULL
  model <- NULL
  exposure <- NULL
  natal <- NULL
  Preterm <- NULL
  
  for (regression in neuroimaging_regressions) {
    beta <- c(beta, regression$coefficients[8,1])
    se <- c(se, regression$coefficients[8,2])
    pvals <- c(pvals, regression$coefficients[8,4])
    r2 <- c(r2, regression$r.squared)
    #exposure <- c()
  }
  
  return(data.frame(beta = beta,
                    se = se,
                    pvals = pvals,
                    r2 = r2,
                    model = "HCA",
                    exposure = "perinatal",
                    natal = "perinatal",
                    Preterm = "Preterm",
                    brain_var = unlist(neuroimaging_list)
                    
  ))
}


#Here tell the function which regressions to use
H9_HCA <- summary_table_h2_preterm(neuroimaging_regressions, neuroimaging_list, neuroimaging_list)
H9_HCA


### H10
### Model 10: BPD
i <- 1

for (brain_var in neuroimaging_list) { 
  linear_mod_output <- summary(lm(scale(brain_volumes_preterms[brain_var]) ~ scale(gest_age) 
                                  + infant_sex 
                                  + scale(BW_Z) 
                                  + scale(gest_scan) 
                                  + Phase
                                  + BPD
                                  + scale(scores), 
                                  data = brain_volumes_preterms))
  rownames(linear_mod_output$coefficients)[8] <- brain_var
  neuroimaging_regressions[[i]] <- linear_mod_output
  i <- (i + 1)
}

########

#MODEL 1 (term : age  + infant_sex) This function outputs my various regression summaries into a table
summary_table_h2_preterm <- function(neuroimaging_regressions, neuroimaging_list, brain_var) {
  
  beta <- NULL
  se <- NULL
  pvals <- NULL
  r2 <- NULL
  model <- NULL
  exposure <- NULL
  natal <- NULL
  Preterm <- NULL
  
  for (regression in neuroimaging_regressions) {
    beta <- c(beta, regression$coefficients[8,1])
    se <- c(se, regression$coefficients[8,2])
    pvals <- c(pvals, regression$coefficients[8,4])
    r2 <- c(r2, regression$r.squared)
    #exposure <- c()
  }
  
  return(data.frame(beta = beta,
                    se = se,
                    pvals = pvals,
                    r2 = r2,
                    model = "BPD",
                    exposure = "postnatal",
                    natal = "postnatal",
                    Preterm = "Preterm",
                    brain_var = unlist(neuroimaging_list)
                    
  ))
}


#Here tell the function which regressions to use
H10_BPD <- summary_table_h2_preterm(neuroimaging_regressions, neuroimaging_list, neuroimaging_list)
H10_BPD

### H11
### Model 11: Sepsis
i <- 1

for (brain_var in neuroimaging_list) { 
  linear_mod_output <- summary(lm(scale(brain_volumes_preterms[brain_var]) ~ scale(gest_age) 
                                  + infant_sex 
                                  + scale(BW_Z) 
                                  + scale(gest_scan) 
                                  + Phase
                                  + Sepsis
                                  + scale(scores), 
                                  data = brain_volumes_preterms))
  rownames(linear_mod_output$coefficients)[8] <- brain_var
  neuroimaging_regressions[[i]] <- linear_mod_output
  i <- (i + 1)
}

########

#MODEL 1 (term : age  + infant_sex) This function outputs my various regression summaries into a table
summary_table_h2_preterm <- function(neuroimaging_regressions, neuroimaging_list, brain_var) {
  
  beta <- NULL
  se <- NULL
  pvals <- NULL
  r2 <- NULL
  model <- NULL
  exposure <- NULL
  natal <- NULL
  Preterm <- NULL
  
  for (regression in neuroimaging_regressions) {
    beta <- c(beta, regression$coefficients[8,1])
    se <- c(se, regression$coefficients[8,2])
    pvals <- c(pvals, regression$coefficients[8,4])
    r2 <- c(r2, regression$r.squared)
    #exposure <- c()
  }
  
  return(data.frame(beta = beta,
                    se = se,
                    pvals = pvals,
                    r2 = r2,
                    model = "Sepsis",
                    exposure = "postnatal",
                    natal = "postnatal",
                    Preterm = "Preterm",
                    brain_var = unlist(neuroimaging_list)
                    
  ))
}


#Here tell the function which regressions to use
H11_Sepsis <- summary_table_h2_preterm(neuroimaging_regressions, neuroimaging_list, neuroimaging_list)
H11_Sepsis

### H12
### Model 12: NEC
i <- 1

for (brain_var in neuroimaging_list) { 
  linear_mod_output <- summary(lm(scale(brain_volumes_preterms[brain_var]) ~ scale(gest_age) 
                                  + infant_sex 
                                  + scale(BW_Z) 
                                  + scale(gest_scan) 
                                  + Phase
                                  + NEC
                                  + scale(scores), 
                                  data = brain_volumes_preterms))
  rownames(linear_mod_output$coefficients)[8] <- brain_var
  neuroimaging_regressions[[i]] <- linear_mod_output
  i <- (i + 1)
}

########

#MODEL 1 (term : age  + infant_sex) This function outputs my various regression summaries into a table
summary_table_h2_preterm <- function(neuroimaging_regressions, neuroimaging_list, brain_var) {
  
  beta <- NULL
  se <- NULL
  pvals <- NULL
  r2 <- NULL
  model <- NULL
  exposure <- NULL
  natal <- NULL
  Preterm <- NULL
  
  for (regression in neuroimaging_regressions) {
    beta <- c(beta, regression$coefficients[8,1])
    se <- c(se, regression$coefficients[8,2])
    pvals <- c(pvals, regression$coefficients[8,4])
    r2 <- c(r2, regression$r.squared)
    #exposure <- c()
  }
  
  return(data.frame(beta = beta,
                    se = se,
                    pvals = pvals,
                    r2 = r2,
                    model = "NEC",
                    exposure = "postnatal",
                    natal = "postnatal",
                    Preterm = "Preterm",
                    brain_var = unlist(neuroimaging_list)
                    
  ))
}


#Here tell the function which regressions to use
H12_NEC <- summary_table_h2_preterm(neuroimaging_regressions, neuroimaging_list, neuroimaging_list)
H12_NEC

### H10
### Model 10: ROP
i <- 1

for (brain_var in neuroimaging_list) { 
  linear_mod_output <- summary(lm(scale(brain_volumes_preterms[brain_var]) ~ scale(gest_age) 
                                  + infant_sex 
                                  + scale(BW_Z) 
                                  + scale(gest_scan) 
                                  + Phase
                                  + ROP
                                  + scale(scores), 
                                  data = brain_volumes_preterms))
  rownames(linear_mod_output$coefficients)[8] <- brain_var
  neuroimaging_regressions[[i]] <- linear_mod_output
  i <- (i + 1)
}

########

#MODEL 1 (term : age  + infant_sex) This function outputs my various regression summaries into a table
summary_table_h2_preterm <- function(neuroimaging_regressions, neuroimaging_list, brain_var) {
  
  beta <- NULL
  se <- NULL
  pvals <- NULL
  r2 <- NULL
  model <- NULL
  exposure <- NULL
  natal <- NULL
  Preterm <- NULL
  
  for (regression in neuroimaging_regressions) {
    beta <- c(beta, regression$coefficients[8,1])
    se <- c(se, regression$coefficients[8,2])
    pvals <- c(pvals, regression$coefficients[8,4])
    r2 <- c(r2, regression$r.squared)
    #exposure <- c()
  }
  
  return(data.frame(beta = beta,
                    se = se,
                    pvals = pvals,
                    r2 = r2,
                    model = "ROP",
                    exposure = "postnatal",
                    natal = "postnatal",
                    Preterm = "Preterm",
                    brain_var = unlist(neuroimaging_list)
                    
  ))
}


#Here tell the function which regressions to use
H13_ROP <- summary_table_h2_preterm(neuroimaging_regressions, neuroimaging_list, neuroimaging_list)
H13_ROP

# Interaction 

### H1 
### Model 1: interaction / fully adjusted
i <- 1
for (brain_var in neuroimaging_list) { 
  linear_mod_output <- summary(lm(scale(brain_volumes_preterms[brain_var]) ~ 
                                    scale(gest_age) 
                                  + infant_sex 
                                  + scale(BW_Z) 
                                  + scale(gest_scan) 
                                  + Phase
                                  # + Mat_Age
                                  + Smoked_in_pregnancy 
                                  # + Ever_smoker 
                                  # + GestDiabetes 
                                  # + steroids
                                  # + MgSO4
                                  + preeclampsia
                                  + HCA
                                  + ROP 
                                  + NEC 
                                  + Sepsis
                                  + BPD
                                  + scale(scores), 
                                  data = brain_volumes_preterms))
  rownames(linear_mod_output$coefficients)[14] <- brain_var
  neuroimaging_regressions[[i]] <- linear_mod_output
  i <- (i + 1)
}

########

#MODEL 1 (term : age  + infant_sex) This function outputs my various regression summaries into a table
summary_table_h1_preterm <- function(neuroimaging_regressions, neuroimaging_list, brain_var) {
  
  beta <- NULL
  se <- NULL
  pvals <- NULL
  r2 <- NULL
  model <- NULL
  exposure <- NULL
  natal <- NULL
  Preterm <- NULL
  
  for (regression in neuroimaging_regressions) {
    beta <- c(beta, regression$coefficients[14,1])
    se <- c(se, regression$coefficients[14,2])
    pvals <- c(pvals, regression$coefficients[14,4])
    r2 <- c(r2, regression$r.squared)
    #exposure <- c()
  }
  
  return(data.frame(beta = beta,
                    se = se,
                    pvals = pvals,
                    r2 = r2,
                    model = "fully adjusted model",
                    exposure = "fully_adj",
                    natal = "postnatal",
                    Preterm = "both",
                    brain_var = unlist(neuroimaging_list)
                    
  ))
}


#Here tell the function which regressions to use
H1_interaction <- summary_table_h1_preterm(neuroimaging_regressions, neuroimaging_list, neuroimaging_list)
H1_interaction


#####
#### TERMS
### Model_Terms_h1
i <- 1
for (brain_var in neuroimaging_list) { 
  linear_mod_output <- summary(lm(scale(brain_volumes_terms[brain_var]) ~ 
                                    scale(gest_age) 
                                  + infant_sex 
                                  + scale(BW_Z) 
                                  + scale(gest_scan) 
                                  + Phase
                                  + scale(scores),
                                  
                                  data = brain_volumes_terms))
  rownames(linear_mod_output$coefficients)[7] <- brain_var
  neuroimaging_regressions[[i]] <- linear_mod_output
  i <- (i + 1)
}



#######
#MODEL 1 (term : age  + infant_sex) This function outputs my various regression summaries into a table
summary_table_Terms_h1 <- function(neuroimaging_regressions, neuroimaging_list, brain_var) {
  
  beta <- NULL
  se <- NULL
  # Upper_CI <- NULL
  # Lower_CI <- NULL
  pvals <- NULL
  r2 <- NULL
  # add_r2 <- NULL
  model <- NULL
  
  for (regression in neuroimaging_regressions) {
    beta <- c(beta, regression$coefficients[7,1])
    se <- c(se, regression$coefficients[7,2])
    # Upper_CI <- (beta + (1.96*se))
    # Lower_CI <- (beta - (1.96*se))
    pvals <- c(pvals, regression$coefficients[7,4])
    r2 <- c(r2, regression$r.squared)
    #  add_r2 <- c(add_r2, regression$r.squared - )
  }
  
  return(data.frame(brain_var = unlist(neuroimaging_list),
                    beta = beta,
                    se = se,
                    #  Upper_CI = Upper_CI,
                    #  Lower_CI = Lower_CI,
                    pvals = pvals,
                    r2 = r2,
                    #   add_r2 = add_r2,
                    model = "Terms [h1]",
                    exposure = "h1",
                    natal = "postnatal",
                    Preterm = "term"
                    
                    
  ))
}

#Here tell the function which regressions to use
Model_Terms_h1 <- summary_table_Terms_h1(neuroimaging_regressions, neuroimaging_list, neuroimaging_list)
Model_Terms_h1




### H14
### Model 14: Maternal age in pregnancy
i <- 1

for (brain_var in neuroimaging_list) { 
  linear_mod_output <- summary(lm(scale(brain_volumes_preterms[brain_var]) ~ 
                                    scale(gest_age) 
                                  + infant_sex 
                                  + scale(BW_Z) 
                                  + scale(gest_scan) 
                                  + Phase
                                  + scale(Mat_Age)
                                  + scale(scores), 
                                  data = brain_volumes_preterms))
  rownames(linear_mod_output$coefficients)[8] <- brain_var
  neuroimaging_regressions[[i]] <- linear_mod_output
  i <- (i + 1)
}

########

#MODEL 1 (term : age  + infant_sex) This function outputs my various regression summaries into a table
summary_table_h14_mat_age <- function(neuroimaging_regressions, neuroimaging_list, brain_var) {
  
  beta <- NULL
  se <- NULL
  pvals <- NULL
  r2 <- NULL
  model <- NULL
  exposure <- NULL
  natal <- NULL
  Preterm <- NULL
  
  for (regression in neuroimaging_regressions) {
    beta <- c(beta, regression$coefficients[8,1])
    se <- c(se, regression$coefficients[8,2])
    pvals <- c(pvals, regression$coefficients[8,4])
    r2 <- c(r2, regression$r.squared)
    #exposure <- c()
  }
  
  return(data.frame(beta = beta,
                    se = se,
                    pvals = pvals,
                    r2 = r2,
                    model = "Maternal age",
                    exposure = "antenatal",
                    natal = "antenatal",
                    Preterm = "Preterm",
                    brain_var = unlist(neuroimaging_list)
                    
  ))
}


#Here tell the function which regressions to use
H14_mat_age <- summary_table_h14_mat_age(neuroimaging_regressions, neuroimaging_list, neuroimaging_list)
H14_mat_age


#####

neuroimaging <- rbind(
  Model_Terms_h1,
  H1,
  H1_interaction,
  H3_smoked_pregnancy,
  H8_preeclampsia,
  H9_HCA,
  H10_BPD,
  H11_Sepsis,
  H12_NEC,
  H13_ROP,
  H14_mat_age
  
  # Preterm_H15_ALL
)

plot <-neuroimaging

plot <- neuroimaging %>% filter(
  brain_var == "Ventricles"|
  brain_var ==   "CSF"|
  brain_var ==   "Brainstem"|
  brain_var ==   "Cerebellum"|
  brain_var ==   "Cortical_grey_matter"|
  brain_var ==   "Hippocampi_and_Amygdala"|
  brain_var ==   "Deep_grey_matter"|
    brain_var ==   "White_matter")


plot <- plot %>% mutate(model = factor(
  model,
  levels = c(
    "BPD",
    "Sepsis",
    "NEC",
    "ROP",
    "HCA",
    "preeclampsia",
    "Smoked during pregnancy",
    "Maternal age",
    "fully adjusted model",
    "gest_age  + infant_sex + BW_Z + gest_scan",
    "Terms [h1]"
  )
))



plot  <- plot  %>% 
  mutate(brain_var =
           case_when(
             brain_var == "Hippocampi_and_Amygdala" ~ "Hippocampi and Amygdalae",
              brain_var == "White_matter" ~ "White matter",
              brain_var == "Deep_grey_matter" ~ "Deep grey matter",
              brain_var == "Cortical_grey_matter" ~ "Cortical grey matter",
              brain_var == "Cerebellum" ~ "Cerebellum",
              brain_var == "Brainstem" ~ "Brainstem",
              brain_var == "CSF" ~ "CSF",
              brain_var == "Ventricles" ~ "Ventricles",
                     TRUE ~ "misc"
         )) %>% 
mutate(brain_var = factor(
  brain_var,
  levels = c(
    # "PSMD",
    #"PSFA",
    "White matter",
    "Deep grey matter",
    "Hippocampi and Amygdalae",
    "Cortical grey matter",
    "Cerebellum",
    "Brainstem",
    "CSF",
    "Ventricles"
  )))

x <- ggplot(plot,
            aes(
              x = brain_var,
              y = beta,
              colour = model,
              group = model,
              shape = model
            )) +
  
  geom_point(position = position_dodge(width = 0.9),
             size = 2.2,
             stroke = 0.9) +
  
  geom_errorbar(
    aes(ymin = beta - (1.96 * se), ymax = beta + (1.96 * se)),
    position = position_dodge(0.9),
    width = 0.3,
    colour = "darkgrey",
    alpha = 0.9,
    size = 0.8
  ) +
  
  geom_hline(yintercept = 0, linetype = "dotted") +
  
  coord_flip() +
  theme_bw() +
  xlab("neuroimaging phenotype") +
  ylab("Standardised effect size")

x +
  theme_classic() +
  theme(
    
    strip.text = element_text(
      size = 8,
      face = "bold",
      family = "sans",
      colour = "black"),
    
    axis.text.x = element_text(
      vjust = 0.5,
      hjust = 1,
      size = 8,
      family = "sans"),
    
    axis.text.y = element_blank(),
    axis.ticks.y=element_blank(),
    
    axis.title.y = element_text(
      size = 10,
      face = "bold",
      family = "sans"),
    
    axis.title.x = element_text(
      size = 10,
      face = "bold",
      family = "sans")
) +
  scale_shape_manual(values = c(12, 
                                8,
                                # 13,
                                # 14,
                                #5,
                                #12,
                                18,
                                2,
                                23,
                                17,
                                4,
                                 6,
                                9,
                                1,
                                16
  )) +
  

  facet_wrap(vars(brain_var), 
             scales = "free_y", 
             nrow = 2, 
             strip.position = "top") +



# 5.56 x 7.80

scale_colour_manual(values = c("#7D1D67",
                            
                               # "#931E6F",
                               "#A82276",
                                "#BD2879",
                               "#D13279",
                               "#E33E76",
                               "#F44C6F",
                               "#F8676C",
                               "#FC7D6C",
                               "#FF9271",
                               "#FFA578",
                              # "#FFB783",
                               #"#FFC891", 
                            "cornflowerblue"
                              # "indianred2"
                            )) +
  ylim(-0.6, 0.6)



####