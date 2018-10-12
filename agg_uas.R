# solusi soal no 1
cor.mtx1 <- matrix(c(1, 0.81, -0.72, 0.81, 1, -0.61, -0.72, -0.61, 1), nrow = 3, byrow = T)
cor.mtx2 <- matrix(c(1, 0.12, -0.23, 0.12, 1, -0.05, -0.23, -0.05, 1), nrow = 3)
eig1 <- eigen(cor.mtx1)
eig2 <- eigen(cor.mtx2)

# solusi soal no 2
library(expm)
pxx_0.5 <- sqrtm(solve(matrix(c(1, 0.37, .21, .37, 1, .35, .21, .35, 1), nrow = 2)))
pyy_1 <- solve(matrix(c(1, 0.8, 0.8, 1), nrow = 2))
pxy <- matrix(c(.26, .67, .34, .33, .59, .34), nrow = 3)
pyx <- t(pxy)
comb <- pxx_0.5 %*% pxy %*% pyy_1 %*% pyx %*% pxx_0.5
eig.res <- eigen(comb)

# solusi soal 3
setwd('/storage/Code/r/apg')
df <- read.csv('CountryData.csv')
count_na <- colSums(apply(df, 2, is.na))
df <- df[,count_na < 45]
cmp_index <- complete.cases(df[,colnames(df)])
df <- df[cmp_index,]
df <- df[, -c(1,2)]

pcomp <- prcomp(df, center = T, scale. = T)
plot(pcomp, type = 'l')
summary(pcomp)
library(ggbiplot)
ggbiplot(pcomp, labels = country_name)

#library used
#library(ggbiplot)
#library(biotools)

#data load
setwd('/storage/Code/r/apg/uas')
df <- read.csv('CountryData.csv')
head(df)
dim(df)

#preprocess
count_na <- colSums(apply(df, 2, is.na))
df <- df[,count_na < 45]
cmp_index <- complete.cases(df[,colnames(df)])
df <- df[cmp_index,]
df <- df[, -c(1,2)]

library(psych)
KMO(cor(df))
df <- df[,!(colnames(df) %in% c('growth', 'death', 'migr', 'inflation', 'gasExp'))]
eig <- eigen(cov(scale(df)))
eig <- eigen(cov(scale(df)))
sum(eig$values >= 1)
fa5 <- factanal(x = df, factors = 5, rotation = 'varimax', lower = 0.03)
print(fa5$loadings, cutoff = 0.5, sort = T)

df.kl <- df[,c('GDP', 'GDPgrowth', 'labor', 'tax', 'budget', 'exports', 'imports')]
clusters <- kmeans(df.kl, 2, iter.max = 10000, nstart = 1000)
df.kl <- data.frame(name = country_name, df.kl, clust = clusters$cluster)

#two-way manova, BoxM
#basic_function2.R
#library for manova is stats
#summary for multivariate, summary.aov for univariate
#boxM use biotools

#components obtained from cov and cor are differ
#relative importance of components is affected by standardization
#rotation of axes maximize variations
#scree plot help to decide how many components included (elbow)
pcomp <- prcomp(df, center = T, scale. = T)
eig.test <- eigen(cor(scale(df))) #give the same eigen val and vec
summary(pcomp)
head(eig.test$vectors, n = 2)
head(pcomp$rotation, n = 2)
plot(pcomp, type = 'l')
ggbiplot(pcomp, labels = country_name)
