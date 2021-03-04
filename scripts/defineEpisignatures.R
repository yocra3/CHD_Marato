#'#################################################################################
#'#################################################################################
#' Define episignatures for CHD groups
#' Adapt algorithm from PMID: 30929737
#' 1. Divide dataset 3/4 train, 1/4 test
#'   - Sample independently each CHD class
#' 2. Run limma of groups each against the rest
#' 3. Select top 1,000 CpGs with minimum difference > 0.1
#' 4. Train SVM with CpGs selected in any comparison
#' 5. Run SVM in test dataset
#' 6. Repeat 4 times, to include all samples in test
#'#################################################################################
#'#################################################################################

## Load libraries
library(minfi)
library(e1071)
library(tidyr)
library(dplyr)
library(limma)

## Load files
load("results/methylation/finalQC_files/v4/gset.autosomic.Rdata")
resFold <- "results/methylation/Episignatures/"
groups <- 4

## Split the dataset in four groups
set.seed(0)
sampTab <- colData(gset)[, c("Sample_Name", "pathClass")] %>%
	data.frame() %>%
	tibble() %>%
	group_by(pathClass) %>%
	mutate(Set = sample(rep(seq_len(groups), each = ceiling(n()/groups)), n()))

## Define the functions
selectDataset <- function(gset, sampTab, batch){
	
	samps <- subset(sampTab, Set != batch)$Sample_Name
	gset[, samps]
}

getFeatures <- function(gset, group){
	
	gset$testVar <- ifelse(gset$pathClass == group, group, "Others")
	model <- model.matrix(~ testVar + Sex, colData(gset))
	
	lmfit <- lmFit(getBeta(gset), design = model)
	lmFite <- eBayes(lmfit)
	tab <- topTable(lmFite, n = Inf, coef = 2)
	tab.fil <- subset(tab, abs(logFC) > 0.1)
	if (nrow(tab.fil) > 1000){
		tab.fil <- tab.fil[1:1000, ]
	}
	rownames(tab.fil)
}

trainSVM <- function(gset){
	
	## Get features distinguishing the groups
	groups <- unique(gset$pathClass)
	featslist <- lapply(groups, getFeatures, gset = gset)
	feats <- unique(unlist(featslist))
	
	mat <- getBeta(gset[feats, ])
	df <- data.frame(pathClass = factor(gset$pathClass), t(mat))
	model_svm <- svm(pathClass ~ ., df)
	list(svm = model_svm, feats = feats)
}

# Compute SVM
## SVM per subgroup ####
#### 3 pathological groups
episig_subset <- lapply(seq_len(groups), function(i){
	train <- selectDataset(gset, sampTab, i)
	test <- gset[, !colnames(gset) %in% colnames(train)]
	
	svm_res <- trainSVM(train)
	
	pred_train <- predict(svm_res$svm, t(getBeta(train[svm_res$feats, ])))
	pred_test <- predict(svm_res$svm, t(getBeta(test[svm_res$feats, ])))
	
	list(svm = svm_res$svm, feats = svm_res$feats, 
		 acc_train = mean(pred_train == train$pathClass),
		 pred_train = pred_train,
		 acc_test = mean(pred_test == test$pathClass),
		 pred_test = pred_test)
})
preds <- unlist(lapply(episig_subset, function(x) x$pred_test))

#### Merge conotruncal malformations and LHLS
episig_subset_mergeClass <- lapply(seq_len(groups), function(i){
	
	gsetmod <- gset
	gsetmod$pathClass <- ifelse(gsetmod$pathClass %in% c("Conotruncal Malformations", "Left heart hypoplasia"),
								"Common malformations", gsetmod$pathClass)
	
	train <- selectDataset(gsetmod, sampTab, i)
	test <- gsetmod[, !colnames(gsetmod) %in% colnames(train)]
	
	svm_res <- trainSVM(train)
	
	pred_train <- predict(svm_res$svm, t(getBeta(train[svm_res$feats, ])))
	pred_test <- predict(svm_res$svm, t(getBeta(test[svm_res$feats, ])))
	
	list(svm = svm_res$svm, feats = svm_res$feats, 
		 acc_train = mean(pred_train == train$pathClass),
		 pred_train = pred_train,
		 acc_test = mean(pred_test == test$pathClass),
		 pred_test = pred_test)
})

## SVM with all samples ####
### 3 CHD groups
svm_all <- trainSVM(gset)
pred_all <- predict(svm_all$svm, t(getBeta(gset[svm_all$feats, ])))

gsetmod <- gset
gsetmod$pathClass <- ifelse(gsetmod$pathClass %in% c("Conotruncal Malformations", "Left heart hypoplasia"),
							"Common malformations", gsetmod$pathClass)
svm_comb <- trainSVM(gsetmod)
save(svm_all, svm_comb, episig_subset_mergeClass, episig_subset, 
	 file = paste0(resFold, "cohort.training.Rdata"))
