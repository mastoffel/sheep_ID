library(caret)

inTrain <- createDataPartition(roh_mod$prop_ROH, p = .10, list = FALSE) %>% 
                as.numeric()
training <- roh_mod[inTrain,]
testing <-  roh_mod[-inTrain,]

mod <- lmer(prop_ROH ~ r_sum_std + het_mean_std + (1|CHR), data = roh_mod)
#mod <- lmer(prop_ROH ~ 1 + r_sum_std + + het_mean_std + (1|CHR), data = roh_mod)

roh_prob <- predict(mod, testing, type = "response")
(cor(roh_prob, testing$prop_ROH))^2
plot(roh_prob, testing$prop_ROH)
