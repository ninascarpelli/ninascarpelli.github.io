#Loading the packages

library(tidyverse)
library(randomForest)

#Setting seed so results can be reproduced

set.seed(123)

#Removing old objects

rm(list = ls())

#Function

getDataPath <- function (...) {
  return(file.path("C:",  ...))
}

#More 

folder <- "AIndices"
step0 <- "9_MergedDF"
step1 <- "10_RF"

#Round1----
# one RF per point
#First thing I'll clean up the file - removing some repeated and unnecessary columns + add the information on the site so next time I run I can just merge the files and do it. I'll create folders named first round; second round; third round inside a RF folder; I'll also add a spreadsheet with the RF results for each round so I can evaluate it later. This part doesn't need to repeat, round 2 should start from line 42 - some bits will be kind of repeated though

# round <- "Round1"
# 
# list_files <- list.files(getDataPath(folder, step0), pattern = "chp", full.names = T)
# 
# file <- list_files[1] 
# 
# labelled <- read.csv(file)
# rownames(labelled) <- labelled$id
# 
# labelled <- select(labelled, everything(), -c(X, component, new_ResultMinute, FileName_start_hyper, FileName_end, files, id_path_hyper)) %>% 
#   separate(id, into = c("site_1", "point_1", "month", "index", "number", "what"), sep = "_", remove = F) %>% 
#   select(distance, reference, fid_what, site, date_time, everything(), -c(site_1, point_1, number, what))
# 
# labelled$index <- as.factor(labelled$index)
# labelled$month <- as.factor(labelled$month)
# labelled$time <- as.integer(labelled$time)
# 
# 
# write.csv(labelled, getDataPath(folder, step1, basename(file)))

#Round2----

round <- "Round2"

list_files <- list.files(getDataPath(folder, step1), pattern = "chp", full.names = T)

#This is about making sure there is the same levels of class labels in the train set as the model data trying to increase accuracy


file <- list_files[10]

labelled <- read.csv(file) %>% 
  select(everything(), -X)
rownames(labelled) <- labelled$id
  
model_data <- filter(labelled, class != "NA") %>%
  select(., class, time, index, month, 22:ncol(labelled)) %>% 
  droplevels(.)


train_index <- sample(1:nrow(model_data), 0.6 * nrow(model_data))
train <- model_data[train_index,]%>%
  droplevels(.)


test_index <- setdiff(1:nrow(model_data), train_index)
test <- model_data[test_index,]%>%
  droplevels(.)




# ct <- ctree(bio ~ ., data=bio, controls = ctree_control(minsplit=30, minbucket=40, maxdepth=5))
# 
# pClassId <- predict(ct)
# 
# # check predicted classes against original class labels
# 
# table(bio$bio, pClassId)
# 
# #accuracy
# 
# (sum(bio$bio==pClassId)) / nrow(bio)
# 
# 
# 
# plot(ct, ip_args=list(pval=FALSE), ep_args=list(digits=0))

set.seed(123)

rf <- randomForest(class ~ ., data = train, importance = T, proximity = T)

print(rf)

varImpPlot(rf)

prediction <- predict(rf, newdata = test)

table(test$class, prediction)


(sum(test$class==prediction))/nrow(test)

#Optimising

importance <- as.data.frame(importance(rf)) %>%
  filter(., MeanDecreaseAccuracy >= 0) %>%
  row.names(.)
# 
model_data <- select(model_data, class, all_of(importance)) %>%
  droplevels(.)

#floor(sqrt(ncol(model_data) - 1))

mtry <- tuneRF(model_data[-1],model_data$class, ntreeTry=500,
               stepFactor=1.5,improve=0.01, trace=TRUE, plot=TRUE)
best.m <- mtry[mtry[, 2] == min(mtry[, 2]), 1]
print(mtry)
print(best.m)

#Optimising

train <- select(train, class, all_of(importance)) %>%
  droplevels(.)

test <- select(test, class, all_of(importance)) %>%
  droplevels(.)

rf <- randomForest(class ~ ., data = train, importance = T, proximity = T, mtry = best.m)

print(rf)

varImpPlot(rf)

prediction <- predict(rf, newdata = test)


table(test$class, prediction)


(sum(test$class==prediction)) / nrow(test)


prediction <- predict(rf, newdata = model_data)


table(model_data$class, prediction)

(sum(model_data$class==prediction)) / nrow(model_data)

#Running it for all the tags

classifier <- select(labelled, class, all_of(importance)) %>% 
  droplevels(.)


# classifier$index <- as.factor(classifier$index)

importance <- as.data.frame(importance(rf))

label_model <- as.data.frame(predict(rf, newdata = classifier))

rownames(label_model) <- labelled$id

label_model <- rename(label_model, "RFclass" = `predict(rf, newdata = classifier)`)

final_df <- cbind(labelled, label_model) %>% 
  select(., class, RFclass, everything())

confusion_matrix <- table(final_df$class, final_df$RFclass)



write.csv(final_df, getDataPath(folder, step1, round, paste(basename(file), "_class_RFlabels.csv", sep = "")), row.names = T)
write.csv(confusion_matrix, getDataPath(folder, step1, round, paste(basename(file), "_ConfusionMatrix.csv", sep = "")), row.names = T)
write.csv(importance, getDataPath(folder, step1, round, paste(basename(file), "importance.csv", sep = "")), row.names = T)
