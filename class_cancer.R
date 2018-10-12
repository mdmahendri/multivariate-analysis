setwd('~/Code/r/apg')
breast_dat <- read.csv('breast_cancer.csv')
head(breast_dat)
table(breast_dat$diagnosis)

#preprocess
breast_dat <- breast_dat[,c(-1,-33)] #remove patient id and blank column
breast_dat$diagnosis <- ifelse(breast_dat$diagnosis == 'M', 1, 0) #M for malignant

#split data
set.seed(123)
train_size <- floor(0.75 * nrow(breast_dat))
train_idx <- sample(seq_len(nrow(breast_dat)), size = train_size)
train_dat <- breast_dat[train_idx,]
test_dat <- breast_dat[-train_idx,]

#first model - using logistic regression
logreg <- glm(
    formula = diagnosis ~ .,
    family = binomial(link = 'logit'),
    data = train_dat,
    maxit = 1000
)
mean(ifelse(logreg$fitted.values >= 0.5, 1, 0) == val_data$diagnosis)

#second model - using knn
#in second model we need to scale features, because it deals with distance
library(class)
sc_train <- scale(train_dat[,-1])
sc_test <- scale(test_dat[,-1])
perf <- numeric(10)
n <- nrow(sc_train)

for (i in 1:length(perf)) {
    kcv <- knn.cv(train = sc_train, cl = train_dat[,1], k = i)
    conf_mtx <- table(kcv, train_dat[,1])
    correct <- sum(diag(conf_mtx))
    perf[i] <- (n-correct) / n
}
which.min(perf) #k = 8
knn8 <- knn(train = sc_train, test = sc_test, cl = train_dat[,1], k = 8)
test_conf <- table(knn8, test_dat[,1])
(nrow(test_dat)-sum(diag(test_conf)))/nrow(test_dat)

#third model - using random forest
#set number of tree as 500
library(randomForest)
randomf <- randomForest(diagnosis ~ ., data = train_dat)

#fourth model - using ann
library(neuralnet)
nn <- neuralnet(
    formula = diagnosis ~ .,
    data = train_dat,
    hidden = length(names(train_dat)),
    algorithm = 'backprop',
    linear.output = F
)

#fifth model - using svm rbf
#do not forget to scale the features
library(e1071)
svm_rbf <- svm(
    formula = diagnosis ~ .,
    data = data.frame(train_dat[,1], sc_train),
    scale = F,
    type = 'C-classification'
)