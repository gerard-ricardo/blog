# load required packages
library(glmmTMB) # for fitting GLMMs
library(ggplot2) # for creating plots
library(dplyr)   # for data manipulation
library(VGAM)    # for transforming data to the logit scale


#1) Loading data
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

table(data1$tot) #display table of tot column

#3) Visualisation
library(ggplot2) #load ggplot2 package
p0 = ggplot()+geom_point(data1, mapping = aes(x = raw.x, y = prop),
                         position = position_jitter(width = .02), alpha = 0.50,
                         size = data1$tot * 0.2) + theme_light() #create scatter plot with jittered points
p0 = p0 + facet_wrap(~factor) + scale_x_log10(name ='raw.x') #add facet and logarithmic x-axis
p0 #display plot

#4) Fitting a generalized linear mixed-effects model
library(lme4) #load lme4 package
library(splines) #load splines package
md3 <- glmer(cbind(suc,(tot - suc)) ~ scale(raw.x) * factor + (1|obs) ,
             family = binomial, data = data1) #fit a generalized linear mixed-effects model
summary(md3) #display model summary

library(RVAideMemoire) #load RVAideMemoire package for overdispersion test
overdisp.glmer(md3) #display overdispersion test result


# 5) fit a GLMM to binary response data
md1 <- glmmTMB(cbind(suc, (tot-suc)) ~ raw.x*factor + (1|obs), data1, family = 'binomial')
# the response is cbind(suc, (tot-suc)), where 'suc' is the number of successes and 'tot' is the total number of trials for each observation
# the predictors are 'raw.x' (a continuous variable) and 'factor' (a categorical variable with three levels)
# '(1|obs)' specifies that each observation has its own random intercept
# 'family = binomial' specifies that the response follows a binomial distribution
summary(md1)
vec.x <- seq(min(data1$raw.x), max(data1$raw.x), length = 100) # create a sequence of values for 'raw.x' to use for plotting
df.x <- expand.grid(raw.x = vec.x, factor = levels(data1$factor))
mm <- model.matrix(~raw.x*factor, df.x)
eta <- mm %*% fixef(md1)$cond # calculate the linear predictor ('eta') on the logit scale for each observation in the data frame
df.x$prediction <- as.vector(exp(eta) / (1 + exp(eta))) # transform the linear predictor back to the response scale to get the predicted probability
se <- sqrt(diag(mm %*% vcov(md1)$cond %*% t(mm))) # calculate the standard error for the prediction at each value of 'raw.x'
# calculate the upper and lower 95% confidence intervals for the prediction
df.x$upper <- exp(eta + 1.96 * se) / (1 + exp(eta + 1.96 * se))
df.x$lower <- exp(eta - 1.96 * se) / (1 + exp(eta - 1.96 * se))
data1$factor = factor(data1$factor, levels = c('a', 'b', 'c')) # reorder the factor levels in the data

#plot the fit
library(ggplot2)
p0= ggplot()
p0= p0+ geom_point(data = data1, aes(x = raw.x, y = prop,alpha = 0.1), color = 'steelblue', size = data1$tot*0.1, position=position_jitter(width = .01))
p0= p0+ geom_line(data = df.x, aes(x = raw.x, y = prediction), color = 'grey30', size=1)
p0= p0+ geom_ribbon(data = df.x, aes(x = raw.x, ymin=lower, ymax=upper,fill='grey'), alpha=0.2)
p0= p0+ scale_x_log10(limits = c(0.9, 100))
p0 = p0+ labs(x=expression(Treatment),
              y=expression(Survival~(prop.)))
p0= p0+ scale_y_continuous( limits = c(-0.05, 1.01))
p0= p0+ theme_light()
p0= p0+ scale_fill_manual( values = c('grey','khaki2'))
p0= p0+ theme(legend.position='none')
p0= p0+ facet_wrap(~factor, nrow = 1)
p0

#6) Getting ECx's
ec = 50 #put in your ecx here
library(dplyr)
group.fac <-df.x %>% group_by(factor)%>%summarise(estimate = max(prediction))%>%as.data.frame()
top1 = group.fac$estimate[1] #the modelled control of factor 1
inhib1 = top1 -((ec/100) * top1) #an x% decrease from the control for factor 1
library(VGAM)
eta1 <- logitlink(inhib1)
data1$factor <- relevel(data1$factor, ref = 'b') #set reference levels for GLMs
md2 <- glmmTMB(cbind(suc,(tot-suc))~raw.x*factor+(1|obs) ,data1,family='binomial' )
data1$factor <- relevel(data1$factor, ref = 'c') #set reference levels for GLMs
md3 <- glmmTMB(cbind(suc,(tot-suc))~raw.x*factor+(1|obs) ,data1,family='binomial' )
data1$factor = factor(data1$factor,levels = c('a', 'b', 'c')) #Set levels in order for model
betas1 = fixef(md1)$cond[1:2] #intercept and slope for ref 1
betas2 = fixef(md2)$cond[1:2] #intercept and slope for ref 2
betas3 = fixef(md3)$cond[1:2]
ecx1 <- (eta1 - betas1[1])/betas1[2]
ecx2 <- (eta1 - betas2[1])/betas2[2]
ecx3 <- (eta1 - betas3[1])/betas3[2]
pd1 <- -cbind(1, ecx1)/betas1[2]
pd2 <- -cbind(1, ecx2)/betas2[2]
pd3 <- -cbind(1, ecx3)/betas3[2]
ff1 = as.matrix(vcov(md1)$cond[1:2,1:2])
ff2 = as.matrix(vcov(md2)$cond[1:2,1:2])
ff3 = as.matrix(vcov(md3)$cond[1:2,1:2])
ec.se1 <- sqrt(((pd1 %*% ff1 )* pd1) %*% c(1, 1))
ec.se2 <- sqrt(((pd2 %*% ff2 )* pd2) %*% c(1, 1))
ec.se3 <- sqrt(((pd3 %*% ff3 )* pd3) %*% c(1, 1))
upper1 = (ecx1+ec.se1*1.96)
lower1 = (ecx1-ec.se1*1.96)
upper2 = (ecx2+ec.se2*1.96)
lower2 = (ecx2-ec.se2*1.96)
upper3 = (ecx3+ec.se2*1.96)
lower3 = (ecx3-ec.se2*1.96)
ec.df1 = data.frame(ecx1, lower1, upper1)
ec.df2 = data.frame(ecx2, lower2, upper2)
ec.df3 = data.frame(ecx3, lower3, upper3)
ecall = cbind(ec.df1, ec.df2, ec.df3)
ec.df1 #this is your factor 1 ecx values
ec.df2 #this is your factor 1 ecx values
ec.df3 #
p0= p0+ geom_hline(yintercept = inhib1, col = 'red', linetype='dashed')
p0


#Visualising the ECX
ecx.all = bind_cols(data.frame(factor = c('a', 'b', 'c')), data.frame(ecx = c(ecx1, ecx2, ecx3)), data.frame(inhib = c(inhib1, inhib1, inhib1)))
upper.all = bind_cols(data.frame(factor = c('a', 'b', 'c')), data.frame(upper = c(upper1, upper2, upper3)), data.frame(inhib = c(inhib1, inhib1, inhib1)))
lower.all = bind_cols(data.frame(factor = c('a', 'b', 'c')), data.frame(lower = c(lower1, lower2, lower3)), data.frame(inhib = c(inhib1, inhib1, inhib1)))
p0 = p0 + geom_point(data = upper.all, aes(x = upper, y = inhib), color = 'red', size=2)
p0 = p0 + geom_point(data = ecx.all, aes(x = ecx, y = inhib), color = 'red', size=2)
p0 = p0 + geom_point(data = lower.all, aes(x = lower, y = inhib), color = 'red', size=2)
p0




