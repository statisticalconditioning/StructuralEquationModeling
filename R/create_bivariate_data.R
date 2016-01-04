#**************************************************************************
# Create bivariate regression data ----------------------------------------
# Author: William Murrah
# Description: Creates data in Table 2.1 of:
#                     Kline, R. B. (2016). Principles and Practice of 
#                         Structural Equation Modeling. guilford, NY.
#**************************************************************************

Case <- LETTERS[1:20]
X <- c(16, 14,16,12,18,18,13,16,18,22,18,19,16,16,22,12,20,14,21,17)
W <- c(48,47,45,45,46,46,47,48,49,49,50,51,52,52,50,51,54,53,52,53)
Y <- c(100,92,88,95,98,101,97,98,110,124,102,115,92,102,104,85,118,105,111,122)

bivariate <- data.frame(Case = Case, X = X, W = W, Y = Y)

save(bivariate, file='data/bivariate.Rdata')
write.csv(bivariate, file = 'data/bivariate.csv', row.names = FALSE)
