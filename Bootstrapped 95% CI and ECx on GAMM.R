# Load packages (if not installed Tools->Install Packages)
library(mgcv)
library(lattice)
library(MASS)
library(MuMIn)

data1 <- read.table(file='https://raw.githubusercontent.com/gerard-ricardo/data/master/mockmulticb', 
                    header= TRUE, dec=',', na.strings=c('','.','NA')) #read data from a URL into R
head(data1) #view first few rows of data
options(scipen = 999) # turn off scientific notation

#2) Organising and wrangling
str(data1) #check data type is correct
data1$raw.x <- as.numeric(as.character(data1$raw.x)) #convert raw.x column to numeric
data1$suc <- as.integer(as.character(data1$suc)) #convert suc column to integer
data1$tot <- as.integer(as.character(data1$tot)) #convert tot column to integer
data1$factor <- as.factor(as.character(data1$factor)) #convert factor column to factor
data1$prop <- data1$suc/data1$tot #calculate prop column by dividing suc column by tot column
data1$obs <- factor(formatC(1:nrow(data1), flag='0', width = 3)) #create a unique tank ID
nrow(data1) #display number of rows in data
str(data1) #check data type of each column
data1 = data1[complete.cases(data1), ] #remove rows with missing values
data1 <- subset(data1, factor!= 'c')  #remove third treatment levels for simplicity


############

m2.gamm<-gamm(cbind(suc,(tot - suc)) ~ s(raw.x, fx = F, k = 3),
              random = list(obs=~1),
              family=binomial,
              data=data1,
              method = "REML") #

df1 <- data.frame(raw.x = seq(min(data1$raw.x), max(data1$raw.x), length = 100)) #setting up  new  data frame (df) defining raw.x values to run 
pred.md1 <- predict(m2.gamm, newdata=df1, type = "response", se.fit=T, lpmatrix=TRUE)
df1$pred = pred.md1$fit
plot(data1$raw.x,(data1$suc / data1$tot),main="GAM") #second plot
lines(df1$raw.x, df1$pred, type = "l", lwd = 2, col = 2, xaxt = "n", las = 1) #plot model mean line
df1$up.ci = pred.md1$fit + 1.96 * pred.md1$se.fit #predicting the data1 95% CI
df1$lo.ci = pred.md1$fit - 1.96 * pred.md1$se.fit #predicting the lower 95% CI
lines(df1$raw.x, df1$up.ci, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot data1 95% CI
lines(df1$raw.x, df1$lo.ci, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot lower 95% CI

####95% CI bootstrap####
set.seed(123)
lst = list() #initial list
fit_gam <- function(data) {
  m2.gamm <- gamm(cbind(suc, (tot - suc)) ~ s(raw.x, fx = F, k = 3),
                  random = list(obs = ~1),
                  family = binomial,
                  data = data,
                  method = "REML")
  pred.md1 <- predict(m2.gamm, newdata=df1, type = "response",se.fit=T)
  out = unname(pred.md1$fit)
  return(out)
} #gam function returns prediction
#bootstrapping
sims = 200
for(i in 1:sims) {
  resampled  = data1[sample(nrow(data1), replace = TRUE), ]
  out = fit_gam(resampled)
  lst[[i]] = out
}
lst  #retuns all predictions in list
df3 <- do.call(rbind, lst)  #add this to 'loop into a list'
df3
#ordering
eta = 0.5*sims
lowerCI = 0.025*sims
upperCI = 0.975*sims
bb_se1<-apply(df3,2,function(X) X[order(X)]) #apply function
df1$boot.pred = bb_se1[eta,] #find the bottom 2.5%
df1$boot_lo = bb_se1[lowerCI,] #find the bottom 2.5%
df1$boot_up = bb_se1[upperCI,] #find the top 2.5%
#plotting
plot(data1$raw.x,(data1$suc / data1$tot),main="GAM") #second plot
lines(df1$raw.x, df1$boot.pred, type = "l", lwd = 2, col = 2, xaxt = "n", las = 1) #plot model mean line
lines(df1$raw.x, df1$boot_lo, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot data1 95% CI
lines(df1$raw.x, df1$boot_up, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot lower 95% CI

#####bootstrap ECX####
#single resampled
ecx = 0.5
resampled  = data1[sample(nrow(data1), replace = TRUE), ]
m2.gamm <- gamm(cbind(suc, (tot - suc)) ~ s(raw.x, fx = F, k = 3),
                random = list(obs = ~1),
                family = binomial,
                data = resampled,
                method = "REML")
pred<- predict(m2.gamm, newdata=df1, type = "response",se.fit=T)
df2 = data.frame(raw.x = df1$raw.x, pred = pred$fit)
ecx.interp2 <- function(ecx, df2) {
  inhibx <- max(df2$pred) * ecx
  nearest_idx <- which.min(abs(df2$pred - inhibx))
  ecx = df2$raw.x[nearest_idx]
  return(ecx)
}
ecx.interp2(0.5, df2)

##loop
set.seed(123)
lst1 = list() #initial list
fit_gam2 <- function(data, ecx) {
  m2.gamm <- gamm(cbind(suc, (tot - suc)) ~ s(raw.x, fx = F, k = 3),
                  random = list(obs = ~1),
                  family = binomial,
                  data = data,
                  method = "REML")
  pred<- predict(m2.gamm, newdata=df1, type = "response",se.fit=T)
  df2 = data.frame(raw.x = df1$raw.x, pred = pred$fit)
  ecx.interp2 <- function(ecx, df2) {
    inhibx <- max(df2$pred) * ecx
    nearest_idx <- which.min(abs(df2$pred - inhibx))
    ecx = df2$raw.x[nearest_idx]
    return(ecx)
  }
  out1 = ecx.interp2(ecx, df2)
  out = unname(out1)
  # lst[[1]] = dd
  # return(lst)
  return(out)
} #gam function returns prediction
fit_gam2(data1, 0.5)
sims = 200
for(i in 1:sims) {
  resampled  = data1[sample(nrow(data1), replace = TRUE), ]
  out = fit_gam2(resampled, ecx)
  lst[[i]] = out
}
lst  #returns all predictions in list
df3 = unlist(lst)
df3
hist(df3)
plot(density(df3)) #density plot
df4 = df3[order(df3)]
df1$boot.pred = df4[eta] #find the bottom 2.5%
df1$boot_lo = df4[lowerCI] #find the bottom 2.5%
df1$boot_up = df4[upperCI] #find the top 2.5%


#Plotting
library(ggplot2)
library(tidybayes)
df5 = data.frame(ecx = df3)
p1 = ggplot(df5, aes(x=ecx))+geom_density(aes(fill = 'steelblue4'), alpha=0.3)+
  stat_pointinterval(aes(y = 0.00, x = ecx),.width = c(.66, .95))
theme_light()
p1 = p1+scale_fill_manual( values = c("steelblue4"))+
  scale_color_manual( values = c("steelblue4"))+theme(legend.position="none")#nice
p1 = p1 + scale_y_continuous(name ="Density")
p1 = p1 + coord_cartesian(ylim = c(0.0, 0.15))
p1


