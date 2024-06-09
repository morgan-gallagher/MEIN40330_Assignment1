library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)

#read in clinical data
df_clinical_patient <- read.delim("data_clinical_patient.txt", row.names = 1, 
                                  comment.char = "#")

#subset dataset to female Black and White
df_clinical_patient <- df_clinical_patient |>
  filter(SEX == "Female") |>
  filter(RACE == "White" | RACE == "Black or African American")

#read in RNA data and clean
RNA_raw <- read.delim("data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt", 
                      check.names = FALSE)
RNA_raw[is.na(RNA_raw)] <- 0
RNA_raw <- RNA_raw[RNA_raw$Hugo_Symbol!='',]
RNA_raw <- RNA_raw[!duplicated(RNA_raw$Hugo_Symbol),]
rownames(RNA_raw) <- RNA_raw$Hugo_Symbol
RNA <- as.data.frame(t(RNA_raw[-1:-2]))
rownames(RNA) <- lapply(rownames(RNA), function(x) {str_sub(x, end = -4)})

# Subset clinical data to match patients with RNA data
clin <- df_clinical_patient[rownames(RNA),]
clin <- clin[grepl("^NA", rownames(clin))==F,]

# Subset RNA data to have patients just in clinical data
RNA <- RNA[rownames(clin), ]

#filter the clinical and RNA data without zero times
clin <- clin[clin$OS_MONTHS > 0, ]
RNA <- RNA[rownames(clin), ]


##############################BOTH RACES########################################

#create survival object with black and white patients
surv_obj <- Surv(time = clin$OS_MONTHS, 
                 event = clin$OS_STATUS=="1:DECEASED")

#overall survival for breast cancer for both races
fit <- survfit(surv_obj ~ 1, data = clin)
ggsurvplot(fit, data = clin, xlab = "Months", ylab = "Survival Probability", 
           title = "Survival - Black and White Patients")

#survival by race
fit_race <- survfit(formula = surv_obj ~ RACE, data = clin)
ggsurvplot(fit_race, data = clin, xlab = "Months", 
           ylab = "Survival Probability", title = "Survival Stratified by Race")

#survival curve differences by race
surv_diff_race <- survdiff(formula = surv_obj ~ RACE, data = clin)


############################SPLIT DATA BY RACE##################################

#divide clinical data by race
clin_black <- clin |>
  subset(RACE == "Black or African American")

clin_white <- clin |>
  subset(RACE == "White")

#divide RNA data by race
RNA_black <- RNA[rownames(clin_black), ]

RNA_white <- RNA[rownames(clin_white), ]


#create survival obj for both races
surv_obj_white <- Surv(time = clin_white$OS_MONTHS, 
                 event = clin_white$OS_STATUS=="1:DECEASED")

surv_obj_black <- Surv(time = clin_black$OS_MONTHS, 
                 event = clin_black$OS_STATUS=="1:DECEASED")


###################################PROGENY######################################

library(progeny)

#progeny for white
zscores_white = as.matrix(t(RNA_white))
pathways_white <- progeny(zscores_white, scale=TRUE, organism="Human")  
path_df_white = as.data.frame(pathways_white)

#progeny for black
zscores_black = as.matrix(t(RNA_black))
pathways_black <- progeny(zscores_black, scale=TRUE, organism="Human")
path_df_black = as.data.frame(pathways_black)


##############################GLMNET BY GENE####################################

library("glmpath")
library("glmnet")
library("penalized")

#top 5 most important genes/pathways
ncut = 5

#glm for black
cvfit_black <- cv.glmnet(data.matrix(RNA_black), surv_obj_black, family="cox", 
                         type.measure = "C")
lambda <- cvfit_black$lambda.min


#get results
cfs_black = coef(cvfit_black,s = lambda)
meaning_coefs_black = rownames(cfs_black)[cfs_black[,1]!= 0]
meaning_vals_black = cfs_black[cfs_black[,1]!=0,]

#get top 5 most important 
vals_surv_black = sort(abs(meaning_vals_black),decreasing = TRUE)[1:ncut]

#print
print(paste(names(vals_surv_black),collapse=" + "))

#run coxph
fit.coxph_black <- coxph(surv_obj_black ~ C3P1 + WDR66 + TMSB15A + CYP2F1 + 
                           SNORA15, data = RNA_black, iter.max = 1000)
ggforest(fit.coxph_black, data = RNA_black, 
         main = "Hazard Ratio - Black Patients")



#glm for white
cvfit_white <- cv.glmnet(data.matrix(RNA_white), surv_obj_white, 
                         family="cox",type.measure = "C")
lambda <- cvfit_white$lambda.min

#get results
cfs_white = coef(cvfit_white,s = lambda)
meaning_coefs_white = rownames(cfs_white)[cfs_white[,1]!= 0]
meaning_vals_white = cfs_white[cfs_white[,1]!=0,]

#get top 5 most important
vals_surv_white = sort(abs(meaning_vals_white),decreasing = TRUE)[1:ncut]

#print string to paste into coxph
print(paste(names(vals_surv_white),collapse=" + "))

# run coxph
fit.coxph_white <- coxph(surv_obj_white ~ OR5H1 + VN1R4 + GALP + CGB8 + 
                           `KRTAP12-1`, data = RNA_white)
ggforest(fit.coxph_white, data = RNA_white, 
         main = "Hazard Ratio - White Patients")


###########################GLMNET BY PATHWAY####################################


#glm for black
cvfit_black <- cv.glmnet(data.matrix(path_df_black), surv_obj_black, 
                         family="cox",type.measure = "C")
lambda <- cvfit_black$lambda.min

#get results
cfs_black = coef(cvfit_black,s=lambda)
meaning_coefs_black = rownames(cfs_black)[cfs_black[,1]!= 0]
meaning_vals_black = cfs_black[cfs_black[,1]!=0,]

#sort the pathways from most important to get top 5
vals_surv_black = sort(abs(meaning_vals_black),decreasing = TRUE)[1:ncut]
print(paste(names(vals_surv_black),collapse=" + "))

#run coxph
fit.coxph_black <- coxph(surv_obj_black ~ VEGF + Estrogen + Hypoxia + p53 + 
                           TGFb, data = path_df_black)
ggforest(fit.coxph_black, data = path_df_black, 
         main = "Hazard Ratio - Black Patients")


#glm for white
cvfit_white <- cv.glmnet(data.matrix(path_df_white), surv_obj_white, 
                         family="cox",type.measure = "C")
lambda <- cvfit_white$lambda.min

#get results
cfs_white = coef(cvfit_white,s=lambda)
meaning_coefs_white = rownames(cfs_white)[cfs_white[,1]!= 0]
meaning_vals_white = cfs_white[cfs_white[,1]!=0,]

#sort to get top 5
vals_surv_white = sort(abs(meaning_vals_white),decreasing = TRUE)[1:ncut]
print(paste(names(vals_surv_white),collapse=" + "))

#run coxph
fit.coxph_white <- coxph(surv_obj_white ~ Estrogen + PI3K + p53 + TNFa + VEGF, 
                         data = path_df_white)
ggforest(fit.coxph_white, data = path_df_white, 
         main = "Hazard Ratio - White Patients")

########################KM Estrogen White#######################################


#get sorted pathway data for estrogen
pathway = 'Estrogen'
pathway_data = path_df_white$Estrogen

#unique values from pathway data
uni_path = sort(unique(pathway_data))

#store results
results_path = matrix(1,length(uni_path))

#record sig of assoc between outcome and pathway data 
for (i in 2:(length(uni_path)-1)){ 
  path_i = 1*(pathway_data>uni_path[i])
  logrank = survdiff(surv_obj_white ~ path_i)
  results_path[i] = logrank$pvalue
}

#get cutoff point with smallest p val
min_p_path = which.min(results_path)

#get optimal threshold for deciding groups
opt_thr = uni_path[min_p_path]

#add column for activity level
path_df_white <- path_df_white %>% 
  mutate(Estrogen_group = ifelse(Estrogen >= 1, "high", 
                        ifelse(Estrogen < 1.001*opt_thr,"low","intermediate")))

#calc number of high/inter/low patients
nhigh = sum(pathway_data>1)
ninter = sum((pathway_data<1) & (pathway_data > 1.001*opt_thr))
nlow = sum(pathway_data<1.001*opt_thr)

#fit KM
KM = survfit(surv_obj_white ~ Estrogen_group,data = path_df_white)

#plot Kaplan Meier 
p <- ggsurvplot(KM, data = path_df_white,pval = TRUE,
                xlab = 'Overall survival time, months',
                
                legend.labs=c(
                  paste("High Estrogen activity,\n",nhigh," patients",sep=""),
                  paste("Intermediate Estrogen activity,\n",ninter," patients",sep=""),
                  paste("Low Estrogen activity,\n",nlow," patients",sep="")),
                
                legend.title="",
                
                title = "Estrogen Activity - White Patients"
)

#customize KM plot
ggpar(p, 
      font.main = c(13, "bold"),
      font.x = c(12, "bold"),
      font.y = c(12, "bold"),
      font.caption = c(12, "bold"), 
      font.legend = c(11, "bold"), 
      font.tickslab = c(13, "bold"))

##############################KM Estrogen White#################################


#get sorted pathway data for estrogen
pathway = 'Estrogen'
pathway_data = path_df_black$Estrogen

#unique values from pathway data
uni_path = sort(unique(pathway_data))

#store results
results_path = matrix(1,length(uni_path))

#record sig of assoc between outcome and pathway data 
for (i in 2:(length(uni_path)-1)){
  path_i = 1*(pathway_data>uni_path[i])
  logrank = survdiff(surv_obj_black ~ path_i)
  results_path[i] = logrank$pvalue
}

#get cutoff point with smallest p val
min_p_path = which.min(results_path)

#get optimal threshold for deciding groups
opt_thr = uni_path[min_p_path]

#add column for activity level
path_df_black <- path_df_black %>% 
  mutate(Estrogen_group = ifelse(Estrogen >= 1, "high", 
                        ifelse(Estrogen < 1.001*opt_thr,"low","intermediate")))

#calc number of high/inter/low patients
nhigh = sum(pathway_data>1)
ninter = sum((pathway_data<1) & (pathway_data > 1.001*opt_thr))
nlow = sum(pathway_data<1.001*opt_thr)

#fit KM
KM = survfit(surv_obj_black ~ Estrogen_group,data = path_df_black)

#plot KM
p <- ggsurvplot(KM, data = path_df_black, pval = TRUE,
              xlab = 'Overall survival time, months',
                
              legend.labs=c(paste("High Estrogen activity,\n",nhigh," patients",sep=""),
                            paste("Intermediate Estrogen activity,\n",ninter," patients",sep=""),
                            paste("Low Estrogen activity,\n",nlow," patients",sep="")),
              
              legend.title="", 
              
              title = "Estrogen Activity - Black Patients"
)

#customize plot
ggpar(p, 
      font.main = c(13, "bold"),
      font.x = c(12, "bold"),
      font.y = c(12, "bold"),
      font.caption = c(12, "bold"), 
      font.legend = c(11, "bold"), 
      font.tickslab = c(13, "bold"))

