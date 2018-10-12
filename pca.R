library(ggbiplot)

setwd('/storage/Code/r/apg')
df <- read.csv('CountryData.csv')
head(df)
dim(df)

#Factor analysis is based on correlations, so we need a lot of sample to get pretty
#accurate estimate of correlations
#As a rule of thumb, use sample size five times the number of variable used
#so in the step below, remove NA in row and column to satisfy requirement above
#and get good correlations matrix
count_na <- colSums(apply(df, 2, is.na))
sum(count_na > 45)
df <- df[,count_na < 45]
cmp_index <- complete.cases(df[,colnames(df)])
df <- df[cmp_index,]
country_name <- df$country
df <- df[, -c(1,2)]

pcomp <- prcomp(df, center = T, scale. = T)
plot(pcomp, type = 'l')
summary(pcomp)
library(ggbiplot)
ggbiplot(pcomp, labels = country_name)

library(psych)
#If standardized measurements are used, we replace S 
#by the sample correlation matrix R.
corr_data <- cor(df)
KMO(corr_data)
df <- df[,!(colnames(df) %in% c('growth', 'death', 'migr', 'inflation', 'gasExp'))]
cortest.bartlett(cor(df), n = 188)
eig <- eigen(cov(scale(df)))
sum(eig$values >= 1)
fa <- factanal(df, factors = 5, rotation = "varimax", lower = 0.05)
fa
#1:pop, GDP, labor, exports, imports, elecProd, elecCons, elecCap, mainlines, cell
#netUsers, roadways
#2:elecImp, petroProd, petroImp, netHosts, airports
#3:birth, infant, life, fert, GDPcapita, 
#4:area, gasProd, gasCons
#5: