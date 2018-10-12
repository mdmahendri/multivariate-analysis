#library used
#library(ggbiplot)
#library(biotools)

#data load
setwd('/storage/Code/r/apg')
df <- read.csv('CountryData.csv')
head(df)
dim(df)

#preprocess
count_na <- colSums(apply(df, 2, is.na))
sum(count_na > 45) #show how many columns have NA more than 45
df <- df[,count_na < 45]
cmp_index <- complete.cases(df[,colnames(df)])
df <- df[cmp_index,]
country_name <- df$country
df <- df[, -c(1,2)]

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
eigen(cor(scale(df))) #give the same eigen val and vec
summary(pcomp)
plot(pcomp, type = 'l')
ggbiplot(pcomp, labels = country_name)
