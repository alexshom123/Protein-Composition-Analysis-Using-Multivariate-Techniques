#read data into R
data <- read.csv('Proteins.csv');

#select random sample
set.seed(10499372); sample(2:ncol(data),size=15,replace=FALSE) 

dat <- data.frame(data$FoodGroup, data$Niacin_mg, data$Phosphorus_mg, data$Fat_g,
                  data$Copper_mcg, data$Protein_g, data$Folate_mcg, data$VitB6_mg,
                  data$Thiamin_mg, data$Zinc_mg, data$VitC_mg, data$Magnesium_mg,
                  data$Selenium_mcg, data$VitA_mcg, data$VitB12_mcg, data$Iron_mg)

#correlation between the variables
cor(dat[,2:16])

#Principal Component Analysis
pca.mod <- princomp(dat[,2:16],cor=TRUE); pca.mod

#eigen values and the proportion of variance explained.
var.pc <- pca.mod$sdev^2; var.pc
prop.var.exp <- var.pc/sum(var.pc); prop.var.exp

#Coefficient of PC
pca.mod$loadings

# improves visualization
library(factoextra)

# generate the principal components
res.pca <- prcomp(dat[1:1244,2:16],  scale = TRUE)

# create a scree plot
fviz_eig(res.pca)

# visualise the clustering within individual samples
fviz_pca_ind(res.pca, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = dat$data.FoodGroup, 
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Status") 

# visualise the relationships among the variables
fviz_pca_var(res.pca, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = dat$data.FoodGroup, 
             col.ind = "contrib", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = FALSE,
             legend.title = "Status") 

# visualise the clustering in the individual samples and variables contributions
fviz_pca_biplot(res.pca, geom.ind = "point", pointshape = 21, 
                pointsize = 2, 
                fill.ind = dat$data.FoodGroup, 
                col.ind = "#696969", 
                palette = "jco", 
                addEllipses = TRUE,
                label = "var",
                col.var = "black",
                repel = TRUE,
                legend.title = "Status") 

#load libraries for factor analysis
library(nFactors)
library(psych)

#parallel analyis to determine n factors
ev <- eigen(cor(dat[-1]));ev
ap <- parallel(subject=nrow(dat[-1]),var=ncol(dat[-1]),rep=100,cent=0.05); ap
ns <- nScree(x=ev$values,aparallel=ap$eigen$qevpea);ns
plotnScree(ns)

#FA with ML extraction method
fa.ml <- factanal(dat,factors=5,rotation="varimax"); fa.ml

#Show only if loadings >0.4
print(fa.ml$loadings,cutoff=0.4,sort=TRUE)

#Check assumptions

#linearity can be checked with a matrix scatterplot 
#Matrix scatter plot
pairs(~.,data=dat[,2:16],bg=as.numeric(dat$FoodGroup)+1,pch=15,upper.panel=NULL)

#Test multivariate normality
library(MVN)
hz.test <- mvn(dat$FoodGroup.demo,mvnTest="hz")
hz.test$multivariateNormality
hz.test$univariateNormality

#Chi-square plot for outliers
mvn(dat[2:16],mvnTest="hz",multivariatePlot="qq")
mvn(dat[2:16],mvnTest="hz",multivariateOutlierMethod="quan",showOutliers=TRUE)

#test assumption of multicolinearity (regress each variable onto the rest)
a <- lm(dat$data.Niacin_mg~.,data=dat[-1])
summary(a)

b <- lm(dat$data.Phosphorus_mg~.,data=dat[-1])
summary(b)

c <- lm(dat$data.Fat_g~.,data=dat[-1])
summary(c)

d <- lm(dat$data.Copper_mcg~.,data=dat[-1])
summary(d)

e <- lm(dat$data.Protein_g~.,data=dat[-1])
summary(e)

f <- lm(dat$data.Folate_mcg~.,data=dat[-1])
summary(f)

g <- lm(dat$data.VitB6_mg~.,data=dat[-1])
summary(g)

h <- lm(dat$data.Thiamin_mg~.,data=dat[-1])
summary(h)

i <- lm(dat$data.Zinc_mg~.,data=dat[-1])
summary(i)

j <- lm(dat$data.Zinc_mg~.,data=dat[-1])
summary(j)

k <- lm(dat$data.VitC_mg~.,data=dat[-1])
summary(k)

l <- lm(dat$data.Magnesium_mg~.,data=dat[-1])
summary(l)

m <- lm(dat$data.Selenium_mcg~.,data=dat[-1])
summary(m)

n <- lm(dat$data.VitA_mcg~.,data=dat[-1])
summary(n)

o <- lm(dat$data.VitB12_mcg~.,data=dat[-1])
summary(o)

p <- lm(dat$data.Iron_mg~.,data=dat[-1])
summary(p)















