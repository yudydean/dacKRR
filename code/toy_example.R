source("functions.R")

## generating simulation data
ngrid = 2e3; dim = 2; cp = 0.4; peaks = 2; 
datagen = DataGenFunc(ngrid, dim = dim, cp = cp, peaks = peaks)
my.data = datagen$my.data
x.test = data.frame(x = datagen$x.grid)

## predicting y test 
RDLRT(my.data,x.test=x.test)
