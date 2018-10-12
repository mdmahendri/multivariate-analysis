setwd('~/Code/data-r/apg')
breast_dat <- read.csv('breast_cancer.csv')
breast_dat <- breast_dat[-33]
sapply(breast_dat, function(x) sum(is.na(x)))

#split dataset
set.seed(123)
train_size <- floor(0.75 * nrow(breast_dat))
train_idx <- sample(seq_len(nrow(breast_dat)), size = train_size)
train_dat <- breast_dat[train_idx,]
train_class <- train_dat[[2]]
train_dat <- train_dat[c(-1,-2)]
test_dat <- breast_dat[-train_idx,]
test_class <- test_dat[[2]]
test_dat <- test_dat[c(-1,-2)]
train_dat <- scale(train_dat)
test_dat <- scale(test_dat)

perf <- numeric(10)
library(class)
for (i in 1:length(perf)) {
    kcv <- knn.cv(train = train_dat, cl = train_class, k = i)
    conf_mtx <- table(kcv, train_class)
    n <- nrow(train_dat)
    correct <- sum(diag(conf_mtx))
    perf[i] <- (n-correct) / n
}
which.min(perf) #k = 8
knn8 <- knn(train = train_dat, test = test_dat, cl = train_class, k = 8)
test_conf <- table(knn8, test_class)
(nrow(test_dat)-sum(diag(test_conf)))/nrow(test_dat)
