source("functions.R")

##########################################
## Real Data One: Health screening
## Data source: https://www.kaggle.com/drateendrajha/health-screening-data
##########################################

data = read.csv("Health Screening Data.csv")
data$age = data$AgeinYr
owlab = which(data$BMICat=="Over Weight")
data1 = data[,c(3:5, 7:16)]
search.na = function(x) sum(is.na(x))
na.vec = apply(data1, 1, search.na)

# remove rows with any NAs
data.new = data1[na.vec==0,]
data.new = data.new[,-12]
names(data.new) = c(paste0("x", 1:11), "y")
my.data = data.new[data.new$y<200 & data.new$y>10,]

ap_hi = my.data$x4
ap_lo = my.data$x5

l1 = which(ap_hi < 90)
l2 = which(ap_lo > 120)
l3 = which(ap_hi<=ap_lo)

my.data = my.data[-c(l1, l2, l3),]
x.test= my.data[,-12]
sigma = 5e-4; breaks = c(0, 18, 25, 30, 50, 100, 200); downsample = TRUE
downsamplesize = 2; nnode = 30; basecut = 50; 

RDLRT(my.data,x.test = x.test, downsample = downsample, sigma = sigma, breaks = breaks,
      downsamplesize = downsamplesize, nnode = nnode, basecut = basecut)


##########################################
## Real Data Two: Melbourne housing
## Data source: https://www.kaggle.com/anthonypino/melbourne-housing-market
##########################################

data = read.csv("Melbourne_housing_FULL.csv",head=T)

data1 = data[,c(1, 3:5, 9:15, 18, 19)]

search.na = function(x) sum(is.na(x))

na.vec = apply(data1, 1, search.na)

# remove rows with any NAs
data.new0 = data1[na.vec==0,]
data.new = data.new0[data.new0$Type == "h",]

sizena = which(data.new$Landsize == 0)
data.new = data.new[-sizena,]

lon = data.new$Longtitude
lat = data.new$Lattitude

d = data.frame(lon = lon, lat = lat)
## lake: lat -37.97, long 144.9
d1 = data.frame(lon = 144.9, lat = -37.97)

my.dist = function(x, y) sqrt(sum((x-y)^2))

data.new$Lakedist = apply(d, 1, my.dist, y = d1)

data.new$Distance = as.numeric(as.character(data.new$Distance))

my.data.ori = data.frame(x1 = data.new$Lattitude, x2 = data.new$Longtitude, x3 = data.new$Distance,  y = pp)
my.data = data.frame(x1 = scale(my.data.ori$x1), x2 = scale(my.data.ori$x2), x3 = scale(my.data.ori$x3), y = my.data.ori$y)
x.test = my.data[,-4]

sigma = 10; breaks = c(0, 3000, 5000, 10000, 30000); downsample = TRUE;
downsamplesize = 200; nnode = 100; 

RDLRT(my.data,x.test = x.test, downsample = downsample, sigma = sigma, breaks = breaks,
      downsamplesize = downsamplesize, nnode = nnode, basecut = basecut)


##########################################
## Real Data Three: Earthquake magnitude prediction
## Data source https://www.kaggle.com/datasets/greegtitan/indonesia-earthquake-data
##########################################

data = read.csv('Earthquake_data.csv')
dat = data[,c(3:6)] %>% na.omit
my.data = cbind(apply(dat[,1:(ncol(dat)-1)], 2, scale),y=dat[,ncol(dat)]) %>% data.frame()
names(my.data) = c(paste0("x", 1:(ncol(dat)-1)), "y")
x.test = my.data[,-4]

sigma=1e-3; breaks = c(0,2,3,4,5,6,7,1e6); downsample = TRUE; 
downsamplesize = 1000; nnode=10; method = 'RF'; 

RDLRT(my.data,x.test = x.test, downsample = downsample, sigma = sigma, breaks = breaks,
      downsamplesize = downsamplesize, nnode = nnode, basecut = basecut)








