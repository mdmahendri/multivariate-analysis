#Load data
dt.train <- rbind(iris3[1:25,,1], iris3[1:25,,2], iris3[1:25,,3])
dt.test <- rbind(iris3[26:50,,1], iris3[26:50,,2], iris3[26:50,,3])
clas <- factor(c(rep('s',25), rep('c', 25), rep('v', 25)))

#KNN
library(class)
knn_res <- knn.cv(train = dt.train, cl = clas, k = 3)
table(knn_res, clas)