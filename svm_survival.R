# # run on server
library(lme4)
library(tidyverse)
library(broom.mixed)
#source("theme_clean.R")
library(snpStats)
library(data.table)
library(furrr)
library(caret)

early_survival_top_snps <- read_delim("output/early_survival_top_snps_train08.txt", delim = " ") 

# take fewer loci
survival_svm_full <- early_survival_top_snps %>% 
                        dplyr::select(id:age2_std, contains("roh")) %>% 
                        select(-mum_id) %>% 
                        mutate(birth_year = as.factor(birth_year),
                               sheep_year = as.factor(sheep_year),
                               id = as.factor(id),
                               twin = as.factor(twin),
                               sex = as.factor(sex)) %>% 
                        mutate(survival = ifelse(survival == 1, "yes", "no")) %>% 
                        mutate(survival = factor(survival))

# do some machine learning
set.seed(3322)
survival_svm <- survival_svm_full %>% 
                        mutate(index = 1:nrow(.)) %>% 
                        #select(survival, id, contains("roh")) %>% 
                        #select(-sheep_year) %>% 
                        #select(-birth_year) %>% 
                        filter(is.na(survival)) %>% 
                        sample_frac(1) 

survival_train <- survival_svm %>%   
        group_by(survival) %>% 
        sample_n(0.8)

survival_test <- survival_svm[-survival_train$index, ]

set.seed(3456)
trainIndex <- createDataPartition(survival_svm$survival, p = .8, 
                                  list = FALSE, 
                                  times = 1)

survival_train <- survival_svm[trainIndex, ]
survival_test <- survival_svm[-trainIndex, ]

## 10-fold CV
# fitControl <- trainControl(
#         method = "repeatedcv",
#         number = 10,
#         ## repeated ten times
#         repeats = 10)


# Control params for SVM
ctrl <- trainControl(
        method = "cv", 
        number = 5, 
        classProbs = TRUE,                 
        summaryFunction = twoClassSummary  # also needed for AUC/ROC
)

library(doParallel)
cl <- makePSOCKcluster(8)
registerDoParallel(cl)
# Tune an SVM
set.seed(5628)  # for reproducibility
surv_svm_auc <- train(
        survival ~ ., 
        data = survival_train, 
        method = "svmRadial",               
       # preProcess = c("center", "scale"),  
        metric = "ROC",  # area under ROC curve (AUC)       
        trControl = ctrl,
        tuneLength = 10,
        na.action = "na.omit"
)
stopCluster(cl)

saveRDS(object = surv_svm_auc, file = "output/svm_train.rds")

# surv_svm_auc <- readRDS("output/svm_train.rds")

# surv_svm_auc
# 
# trellis.par.set(caretTheme())
# plot(surv_svm_auc)
# # 
# pred_surv <- predict(surv_svm_auc, survival_test, type = "prob") %>% 
#          mutate(surv = ifelse(no > yes, "no", "yes"))
# # 
# survival_test <- survival_test %>% 
#         drop_na()
# survival_test$pred_surv <- pred_surv$surv
# out <- survival_test %>% select(survival, pred_surv) %>% mutate(pred_surv = factor(pred_surv))
# 
# confusionMatrix(data = out$pred_surv, reference = out$survival)
# head(surv_svm_auc$pred)
