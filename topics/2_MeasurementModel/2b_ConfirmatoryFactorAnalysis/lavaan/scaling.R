gg#**************************************************************************
# Scaling Measurement Models ----------------------------------------------
# Author: William Murrah
# Description: This code explores the methods for scaling measurement 
#              models in SEM. I draw heavily from:
#              Little, P. T. D. (2013). Longitudinal structural equation 
#                modeling. Guilford Press.
#**************************************************************************
# Packages used -----------------------------------------------------------
library(lavaan)

cor <- "
1.00000
0.55226    1.00000
0.56256    0.60307    1.00000
0.31889    0.35898    0.27757    1.00000
0.24363    0.35798    0.31889    0.56014    1.00000
0.32217    0.36385    0.32072    0.56164    0.59738    1.00000"

cormat <- getCov(cor, lower = TRUE, 
                 names = c('glad1', 'cheer1', 'happy1',
                           'glad2', 'cheer2', 'happy2'))


## Marker variable method

model1 <- '
PA =~ glad1 + cheer1 + happy1'
fit1 <- cfa(model1, sample.cov = cormat, sample.nobs = 823)
summary(fit1)

model2 <- '
PA =~ NA*glad1 + 1*cheer1 + happy1'
fit2 <- cfa(model2, sample.cov = cormat, sample.nobs = 823)
summary(fit2)

model3 <- '
PA =~ NA*glad1 + cheer1 + 1*happy1'
fit3 <- cfa(model3, sample.cov = cormat, sample.nobs = 823)
summary(fit3)

model4 <- '
PA =~ NA*glad1 + cheer1 + happy1
PA ~~ 1*PA'
fit4 <- cfa(model4, sample.cov = cormat, sample.nobs = 823)
summary(fit4)

model5 <- '
PA =~ L1*glad1 + L2*cheer1 + L3*happy1
PA ~~ PA
glad1 ~~ glad1
cheer1 ~~ cheer1
happy1 ~~ happy1
L1 == 3 - L2 - L3
'
fit5 <- lavaan(model5, sample.cov = cormat, sample.nobs = 823,
               auto.fix.first=F)
summary(fit5)
