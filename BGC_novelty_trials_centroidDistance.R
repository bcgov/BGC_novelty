# Developing a method for novelty detection in biogeoclimatic projections
# Colin Mahony colin.mahony@gov.bc.ca

library(climr)
library(terra)
library(data.table)
library(bcmaps)
library(ccissr)
library(ranger)
library(scales)
library(EnvStats)
library("FactoMineR")

#BGC model
load("C:/Users/CMAHONY/OneDrive - Government of BC/Data/BGC_models/BGC_RFresp.Rdata") ##load RF model
pred_vars <- BGC_RFresp[["forest"]][["independent.variable.names"]] ##required predictors

# lists
monthcodes <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
elements <- c("Tmin", "Tmax", "PPT")
month.abb.lowercase <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")

bc <- vect(bc_bound())
bc <- project(bc, "EPSG:4326")

# PRISM DEM
# dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/PRISM_dem/", sep="")
dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/PRISM_dem/", sep="")
dem <- rast(paste(dir, "PRISM_dem.asc", sep=""))
dem <- aggregate(dem, fact=3)
dem <- mask(dem, bc)
dem <- trim(dem)
plot(dem)

# climr data
grid <- as.data.frame(dem, cells = TRUE, xy = TRUE)
colnames(grid) <- c("id", "lon", "lat", "elev") # rename column names to what climr expects
clim.grid <- downscale(xyz = grid, 
                       gcms = list_gcms()[1], 
                       ssps = list_ssps()[2],
                       gcm_periods = list_gcm_periods()[3], 
                       run_nm = list_runs_ssp(list_gcms()[1], list_ssps()[2])[3], # not working atm; i've created a github issue
                       # max_run = 2,
                       vars = list_vars()
                       )
#select the run i want (since run.nm isn't working)
clim.grid <- clim.grid[RUN == list_runs_ssp(list_gcms()[1], list_ssps()[2])[3],]
addVars(clim.grid)
clim.grid <- clim.grid[is.finite(CMD.total)] #remove NA rows to have complete cases for RF model

# plot
e=1
m=1
var <- paste(c("Tmin", "Tmax", "PPT")[e], monthcodes[m], sep="_")
X <- rast(dem) # use the DEM as a template raster
X[clim.grid[, id]] <- clim.grid[,..var]
plot(X)

# BGC projections
bgc.pred <- predict(BGC_RFresp, data = clim.grid)[['predictions']]
X <- rast(dem) # use the DEM as a template raster
X[clim.grid[, id]] <- bgc.pred
plot(X)


#historical climate for training points
pts <- fread("C:/Users/CMAHONY/OneDrive - Government of BC/Data/BGC_models/WNA_v13_50-200filtpts_15Nov.csv")
colnames(pts) <- c("id", "BGC", "lon", "lat", "elev") # rename column names to what climr expects
clim.pts <- downscale(xyz = pts,
                      vars = list_vars())
addVars(clim.pts)

# Calculate the centroid climate for the training points
clim.pts.mean <- clim.pts[, lapply(.SD, mean), by = pts$BGC, .SDcols = -c(1,2)]

# check that the centroids make sense
x <- clim.pts.mean$MAT
y <- log(clim.pts.mean$MAP)
plot(x, y, col="white")
text(x,y,clim.pts.mean$pts, cex=0.5)


#---------------------------
# sigma dissimilarity metric methods development
#---------------------------

# example of the Cwhvm1
var1 <- "MAT"
var2 <- "MAP"
bgc.focal <- "CWHmh_OR"
x.fut <- clim.grid[bgc.pred==bgc.focal, ..var1]
y.fut <- clim.grid[bgc.pred==bgc.focal, ..var2]
x <- clim.pts[pts$BGC==bgc.focal, ..var1]
y <- clim.pts[pts$BGC==bgc.focal, ..var2]
x.mean <- clim.pts.mean[clim.pts.mean$pts==bgc.focal, ..var1]
y.mean <- clim.pts.mean[clim.pts.mean$pts==bgc.focal, ..var2]
plot(x,y, xlim=range(c(x, x.fut)), ylim=range(c(y, y.fut)))
points(x.mean, y.mean, pch=16, cex=2)
points(x.fut, y.fut, col="blue", pch=16 )

# z-standardize the data for M distance calculation
bgc.focal <- "SWBvks"
vars <- as.vector(outer(elements, c("wt", "sm"), paste, sep = "_"))
clim.analog <- clim.pts[pts$BGC==bgc.focal, ..vars]
clim.mean <- clim.analog[, lapply(.SD, mean, na.rm = TRUE)]
clim.sd <- clim.analog[, lapply(.SD, sd, na.rm = TRUE)]
clim.target <- clim.grid[bgc.pred==bgc.focal, ..vars]
# z-standardize the analog and target
clim.analog[, (names(clim.analog)) := lapply(names(clim.analog), function(col) {
  (get(col) - unlist(clim.mean)[col]) / unlist(clim.sd)[col]
})]
clim.target[, (names(clim.target)) := lapply(names(clim.target), function(col) {
  (get(col) - unlist(clim.mean)[col]) / unlist(clim.sd)[col]
})]

# M distance in 2 variables. centroid is a vector of zeros due to z-standardization. square root since mahalanobis function returns squared distances
md.analog <- (mahalanobis(clim.analog[,c(4,6)], rep(0, 2), var(clim.analog[,c(4,6)])))^0.5
md.target <- (mahalanobis(clim.target[,c(4,6)], rep(0, 2), var(clim.analog[,c(4,6)])))^0.5
hist(md.target)
hist(md.analog, col=alpha("blue", 0.4), add=T)

# plot to visualize the distances
a <- clim.analog[,c(4,6)]
b <- clim.target[,c(4,6)]
plot(a, xlim=range(c(a[,1], b[,1])), ylim=range(c(a[,2], b[,2])))
points(b, bg=grey(rescale(md.target, to = c(0, 1))), pch=21)

# sigma dissimilarity
sigmaD <- function(md, df){
  p <- pchi(md,df) # percentile of the nearest neighbour distance on the chi distribution with degrees of freedom equaling the dimensionality of the distance measurement (PCs)
  q <- qchi(p,1) # values of the chi percentiles on a standard half-normal distribution (chi distribution with one degree of freedom)
  q[!is.finite(q)] <- max(q[is.finite(q)]) # set infinite values (outside the decimal precision of pchi) to max of finite values
  return(q)
}
sigma.target <- sigmaD(md.target,2)
sigma.analog <- sigmaD(md.analog,2)
hist(sigma.target)
hist(sigma.analog, col=alpha("blue", 0.4), add=T)

# plot to visualize sigma dissimilarity
a <- clim.analog[,c(4,6)]
b <- clim.target[,c(4,6)]
plot(a, col="dodgerblue", xlim=range(c(a[,1], b[,1])), ylim=range(c(a[,2], b[,2])))
points(b, bg=grey(rescale(sigma.target, to = c(0, 1))), pch=21, cex=1.5)

# M distance in all variables. centroid is a vector of zeros due to z-standardization. square root since mahalanobis function returns squared distances
md.analog <- (mahalanobis(clim.analog, rep(0, length(clim.mean)), var(clim.analog)))^0.5
md.target <- (mahalanobis(clim.target, rep(0, length(clim.mean)), var(clim.analog)))^0.5
hist(md.target, xlim=c(0,max(md.target)))
hist(md.analog, col=alpha("blue", 0.4), add=T)

# sigma dissimilarity
sigma.target <- sigmaD(md.target,length(clim.mean)) # percentile of the nearest neighbour distance on the chi distribution with degrees of freedom equaling the dimensionality of the distance measurement (PCs)
sigma.analog <- sigmaD(md.analog,length(clim.mean)) # percentile of the nearest neighbour distance on the chi distribution with degrees of freedom equaling the dimensionality of the distance measurement (PCs)
hist(sigma.target)
hist(sigma.analog, col=alpha("blue", 0.4), add=T)

# plot to visualize sigma dissimilarity
a <- clim.analog[,c(4,6)]
b <- clim.target[,c(4,6)]
plot(a, col="dodgerblue", xlim=range(c(a[,1], b[,1])), ylim=range(c(a[,2], b[,2])))
points(b, bg=grey(rescale(sigma.target, to = c(0, 1))), pch=21, cex=1.5)

# #---------------------------
# # function to calculate sigma dissimilarity from raw data (abandoned because the PCA approach produces the same result with more control)
# #---------------------------
# 
# sigmaD <- function(clim.target, clim.analog){
#   # z-standardize the analog and target
#   clim.mean <- clim.analog[, lapply(.SD, mean, na.rm = TRUE)]
#   clim.sd <- clim.analog[, lapply(.SD, sd, na.rm = TRUE)]
#   clim.analog[, (names(clim.analog)) := lapply(names(clim.analog), function(col) {
#     (get(col) - unlist(clim.mean)[col]) / unlist(clim.sd)[col]
#   })]
#   clim.target[, (names(clim.target)) := lapply(names(clim.target), function(col) {
#     (get(col) - unlist(clim.mean)[col]) / unlist(clim.sd)[col]
#   })]
#   # M distance in all variables. centroid is a vector of zeros due to z-standardization. square root since mahalanobis function returns squared distances
#   md <- (mahalanobis(clim.target, rep(0, length(clim.mean)), var(clim.analog)))^0.5
#   p <- pchi(md,length(clim.mean)) # percentile of the nearest neighbour distance on the chi distribution with degrees of freedom equaling the dimensionality of the distance measurement (PCs)
#   q <- qchi(p,1) # values of the chi percentiles on a standard half-normal distribution (chi distribution with one degree of freedom)
#   q[!is.finite(q)] <- 8 # set infinite values to 8 sigma (outside the decimal precision of pchi) 
#   return(q)
# }
# 
# # test 
# vars <- as.vector(outer(elements, c("wt", "sm"), paste, sep = "_"))
# clim.analog <- clim.pts[pts$BGC==bgc.focal, ..vars]
# clim.target <- clim.grid[bgc.pred==bgc.focal, ..vars]
# sigma.target <- sigmaD(clim.target[,c(4,6)], clim.analog[,c(4,6)])
# a <- clim.analog[,c(4,6)]
# b <- clim.target[,c(4,6)]
# plot(a, col="dodgerblue", xlim=range(c(a[,1], b[,1])), ylim=range(c(a[,2], b[,2])))
# points(b, bg=grey(rescale(sigma.target, to = c(0, 1))), pch=21, cex=1.5)

#---------------------------
# function to calculate sigma dissimilarity that allows truncating minor PCs
#---------------------------
sigmaD <- function(clim.target, clim.analog, pcs = NULL, threshold=.95){

  ## data handling
  clim.analog <- clim.analog[complete.cases(clim.analog)] # remove rows without data
  clim.analog <- clim.analog[, .SD, .SDcols = which(sapply(clim.analog, function(x) var(x, na.rm = TRUE) > 0))]  # Remove zero-variance columns
  clim.target <- clim.target[, .SD, .SDcols = names(clim.analog)]
  
  ## PCA on scaled data
  pca <- prcomp(clim.analog, scale=TRUE)
  pcs.analog <- predict(pca, clim.analog)
  pcs.target <- predict(pca, clim.target)
  pcs.centroid <- apply(pcs.analog, 2, mean)
  
  if(is.null(pcs)){
  ## select number of pcs
  cumvar <- cumsum(pca$sdev^2 / sum(pca$sdev^2)) # vector of cumulative variance explained
    pcs <- which(cumvar >= threshold)[1]
    if(pcs<3) pcs <- 3
  }
  
  md <- (mahalanobis(pcs.target[,1:pcs], pcs.centroid[1:pcs], var(pcs.analog[,1:pcs])))^0.5
  p <- pchi(md,pcs) # percentile of the nearest neighbour distance on the chi distribution with degrees of freedom equaling the dimensionality of the distance measurement (PCs)
  q <- qchi(p,1) # values of the chi percentiles on a standard half-normal distribution (chi distribution with one degree of freedom)
  q[!is.finite(q)] <- 8 # set infinite values to 8 sigma (outside the decimal precision of pchi) 
  
  return(q)
}

# test 
vars <- as.vector(outer(elements, c("wt", "sm"), paste, sep = "_"))
par(mar=c(1,1,1,1), mfrow=c(2,2))
for(i in c(2,3,4,6)){
pcs <- i
clim.analog <- clim.pts[pts$BGC==bgc.focal, ..vars]
clim.target <- clim.grid[bgc.pred==bgc.focal, ..vars]
sigma.target <- sigmaD(clim.target, clim.analog, pcs)
pca <- prcomp(clim.analog, scale=TRUE)
a <- predict(pca, clim.analog)[,1:2]
b <- predict(pca, clim.target)[,1:2]
plot(a, col="dodgerblue", xlim=range(c(a[,1], b[,1])), ylim=range(c(a[,2], b[,2])))
points(b, bg=ColScheme[cut(sigma.target, breakpoints)], pch=21, cex=1.5)
points(a, col="dodgerblue", pch=16)
mtext(paste(bgc.focal, "\n", pcs, "PCs"), line=-2.5, adj = 0.05, )
}

#---------------------------
# scatterplots of all predicted units 
#---------------------------

breakseq <- c(0,4,8)
breakpoints <- c(seq(breakseq[1], breakseq[3], 0.01),199); length(breakpoints)
ColScheme <- c(colorRampPalette(c("gray90", "gray50", "#FFF200", "#CD0000", "black"))(length(breakpoints)))

vars <- as.vector(outer(elements, c("wt", "sm"), paste, sep = "_"))
bgcs.target <- names(table(bgc.pred)[which(table(bgc.pred)>0)]) # list of analogs used in the BGC projection
sigma <- rep(NA, length(bgc.pred)) # initiate a vector to store the sigma dissimilarities
pcs <- 2
par(mar=c(1,1,1,1), mfrow=c(2,2))
for(bgc.focal in bgcs.target){
  clim.analog <- clim.pts[pts$BGC==bgc.focal, ..vars]
  clim.target <- clim.grid[bgc.pred==bgc.focal, ..vars]
  sigma[bgc.pred==bgc.focal] <- sigmaD(clim.target, clim.analog, pcs)
  
  if(dim(clim.target)[1]>25){
    pca <- prcomp(clim.analog, scale=TRUE)
    a <- predict(pca, clim.analog)[,1:2]
    b <- matrix(predict(pca, clim.target)[,1:2], ncol=2)
    plot(a, col="white", xlim=range(c(a[,1], b[,1])), ylim=range(c(a[,2], b[,2])), xaxt="n", yaxt="n")
    points(b, bg=ColScheme[cut(sigma[bgc.pred==bgc.focal], breakpoints)], pch=21, cex=1.5)
    points(a, bg="dodgerblue", pch=21, col="blue", cex=2)
    mtext(bgc.focal, line=-1.5, adj = 0.05, )
  }
  print(paste0(round(which(bgcs.target==bgc.focal)/length(bgcs.target)*100),"%"))
  
}

#---------------------------
# maps at increasing PCs
#---------------------------

vars <- as.vector(outer(elements, c("wt", "sp", "sm", "at"), paste, sep = "_"))
bgcs.target <- names(table(bgc.pred)[which(table(bgc.pred)>0)]) # list of analogs used in the BGC projection
levels(bgc.pred)[-which(levels(bgc.pred)%in%pts$BGC)]
bgcs.target <- bgcs.target[which(bgcs.target%in%pts$BGC)] # BAFAunp not in the training set for some reason
sigma <- rep(NA, length(bgc.pred)) # initiate a vector to store the sigma dissimilarities
par(mar=c(1,1,1,1), mfrow=c(2,2))
for(pcs in c(2,3,4,6)){
  for(bgc.focal in bgcs.target){
    clim.analog <- clim.pts[pts$BGC==bgc.focal, ..vars]
    clim.analog <- clim.analog[complete.cases(clim.analog)]
    clim.analog <- clim.analog[, .SD, .SDcols = which(sapply(clim.analog, function(x) var(x, na.rm = TRUE) > 0))]  # Remove zero-variance columns
    clim.target <- clim.grid[bgc.pred==bgc.focal, .SD, .SDcols = names(clim.analog)]
    sigma[bgc.pred==bgc.focal] <- sigmaD(clim.target, clim.analog, pcs)
  }
  X[clim.grid[, id]] <- sigma
  plot(X, col=ColScheme, axes=F) 
  mtext("Sigma novelty", side=4, line=-4.5, adj = 0.5, font=2)
  mtext(paste0("Base variables", "\n", pcs, "PCs"), line=-3.5, adj = 0.975, )
}

par(mar=c(1,1,1,1), mfrow=c(1,1))
pcs <- 4
for(bgc.focal in bgcs.target){
  clim.analog <- clim.pts[pts$BGC==bgc.focal, ..vars]
  clim.target <- clim.grid[bgc.pred==bgc.focal, ..vars]
  sigma[bgc.pred==bgc.focal] <- sigmaD(clim.target, clim.analog, pcs)
}
X[clim.grid[, id]] <- sigma
plot(X, col=ColScheme) 
mtext(paste(pcs, "PCs"), line=-1.5, adj = 0.95, )


#---------------------------
# scatterplots of all predicted units 
#---------------------------

breakseq <- c(0,4,8)
breakpoints <- c(seq(breakseq[1], breakseq[3], 0.01),199); length(breakpoints)
ColScheme <- c(colorRampPalette(c("gray90", "gray50", "#FFF200", "#CD0000", "black"))(length(breakpoints)))

vars <- pred_vars
bgcs.target <- names(table(bgc.pred)[which(table(bgc.pred)>0)]) # list of analogs used in the BGC projection
sigma <- rep(NA, length(bgc.pred)) # initiate a vector to store the sigma dissimilarities
pcs <- 2
par(mar=c(1,1,1,1), mfrow=c(2,2))
# for(bgc.focal in bgcs.target){
  clim.analog <- clim.pts[pts$BGC==bgc.focal, ..vars]
  clim.analog <- clim.analog[complete.cases(clim.analog)]
  clim.analog <- clim.analog[, .SD, .SDcols = which(sapply(clim.analog, function(x) var(x, na.rm = TRUE) > 0))]  # Remove zero-variance columns
  clim.target <- clim.grid[bgc.pred==bgc.focal, .SD, .SDcols = names(clim.analog)]
  sigma[bgc.pred==bgc.focal] <- sigmaD(clim.target, clim.analog, pcs)
  
  if(dim(clim.target)[1]>50){
    pca <- prcomp(clim.analog, scale=TRUE)
    a <- predict(pca, clim.analog)[,1:2]
    b <- matrix(predict(pca, clim.target)[,1:2], ncol=2)
    plot(a, col="white", xlim=range(c(a[,1], b[,1])), ylim=range(c(a[,2], b[,2])), xaxt="n", yaxt="n")
    points(b, bg=ColScheme[cut(sigma[bgc.pred==bgc.focal], breakpoints)], pch=21, cex=1.5)
    points(a, bg="dodgerblue", pch=21, col="blue", cex=2)
    mtext(paste(bgc.focal, "\n", pcs, "PCs"), line=-2.5, adj = 0.05, )
  }
  print(paste0(round(which(bgcs.target==bgc.focal)/length(bgcs.target)*100),"%"))
  
# }


#---------------------------
# 3D rotatable plot
#---------------------------

bgc.focal <- "CWHxm_WA" # good example of complete novelty
  bgc.focal <- "ESSFwm3" # good example of close analog plus novelty
  bgc.focal <- "CRFdh_CA"
  pcs <- 3

clim.analog <- clim.pts[pts$BGC==bgc.focal, ..vars]
clim.target <- clim.grid[bgc.pred==bgc.focal, .SD, .SDcols = names(clim.analog)]
q <- sigmaD(clim.target, clim.analog, pcs=pcs)

# Perform PCA
pca <- prcomp(clim.analog, scale = TRUE)

# Extract the first three principal components
a <- predict(pca, clim.analog)[, 1:3]
b <- predict(pca, clim.target)[, 1:3]

# Define colors for points in 'b'
b_colors <- ColScheme[cut(q, breakpoints)]

# Create the 3D scatterplot
plot <- plot_ly() %>%
  add_trace(
    x = a[, 1], y = a[, 2], z = a[, 3],
    type = "scatter3d", mode = "markers",
    marker = list(size = 5, color = "dodgerblue", opacity = 1),
    name = "Analog Points"
  ) %>%
  add_trace(
    x = b[, 1], y = b[, 2], z = b[, 3],
    type = "scatter3d", mode = "markers",
    marker = list(size = 6, color = b_colors, opacity = 1),
    name = "Target Points"
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    ),
    title = list(text = bgc.focal, x = 0.05)
  )

# Display the plot
plot

#---------------------------
# maps at increasing PCs
#---------------------------

par(mar=c(1,1,1,1), mfrow=c(2,2))
for(pcs in c(2,3,4,6)){
  for(bgc.focal in bgcs.target){
    clim.analog <- clim.pts[pts$BGC==bgc.focal, ..vars]
    clim.target <- clim.grid[bgc.pred==bgc.focal, ..vars]
    sigma[bgc.pred==bgc.focal] <- sigmaD(clim.target, clim.analog, pcs)
  }
  X[clim.grid[, id]] <- sigma
  plot(X, col=ColScheme, axes=F) 
  mtext("Sigma novelty", side=4, line=-4.5, adj = 0.5, font=2)
  mtext(paste0("RF variables", "\n", pcs, "PCs"), line=-3.5, adj = 0.975, )
}

par(mar=c(1,1,1,1), mfrow=c(1,1))
pcs <- 4
for(bgc.focal in bgcs.target){
  clim.analog <- clim.pts[pts$BGC==bgc.focal, ..vars]
  clim.analog <- clim.analog[complete.cases(clim.analog)]
  clim.analog <- clim.analog[, .SD, .SDcols = which(sapply(clim.analog, function(x) var(x, na.rm = TRUE) > 0))]  # Remove zero-variance columns
  clim.target <- clim.grid[bgc.pred==bgc.focal, .SD, .SDcols = names(clim.analog)]
  sigma[bgc.pred==bgc.focal] <- sigmaD(clim.target, clim.analog, pcs)
}
X[clim.grid[, id]] <- sigma
plot(X, col=ColScheme) 
mtext(paste(pcs, "PCs"), line=-2.5, adj = 0.9, )

#---------------------------
#---------------------------
# alternative function that selects principal components of the combined variation in the analog and target. 
# the principle here is that the metric should capture modes of climate change as well as spatial variation
#---------------------------
#---------------------------

bgc.focal <- "CWHxm_WA" # good example of complete novelty
bgc.focal <- "CRFdh_CA" # good example of complete novelty
vars <- pred_vars

clim.analog <- clim.pts[pts$BGC==bgc.focal, ..vars]
clim.target <- clim.grid[bgc.pred==bgc.focal, ..vars]

sigmaD_spattemp <- function(clim.target, clim.analog, threshold = 0.95, plot2d = FALSE, plot3d = FALSE){
  
  ## data cleaning
  clim.analog <- clim.analog[complete.cases(clim.analog)] # remove rows without data
  clim.analog <- clim.analog[, .SD, .SDcols = which(sapply(clim.analog, function(x) var(x, na.rm = TRUE) > 0))]  # Remove zero-variance columns
  clim.target <- clim.target[, .SD, .SDcols = names(clim.analog)]
  
  ## scale the data to the variance of the analog, since this is what we will ultimately be measuring the M distance in. 
    clim.mean <- clim.analog[, lapply(.SD, mean, na.rm = TRUE)]
    clim.sd <- clim.analog[, lapply(.SD, sd, na.rm = TRUE)]
    clim.analog[, (names(clim.analog)) := lapply(names(clim.analog), function(col) {
      (get(col) - unlist(clim.mean)[col]) / unlist(clim.sd)[col]
    })]
    clim.target[, (names(clim.target)) := lapply(names(clim.target), function(col) {
      (get(col) - unlist(clim.mean)[col]) / unlist(clim.sd)[col]
    })]
  
  ## PCA on pooled target and analog
  s <- sample(1:dim(clim.target)[1], dim(clim.analog)[1], replace = TRUE) # select a random sample of the target population to match the analog points. bootstrap if target population is smaller than analog points
  clim.target.sample <- clim.target[s,]
  pca <- prcomp(rbind(clim.analog, clim.target.sample), scale=FALSE)
  pcs.analog <- data.table(predict(pca, clim.analog))
  pcs.target <- data.table(predict(pca, clim.target))
  
  ## select number of pcs
  cumvar <- cumsum(pca$sdev^2 / sum(pca$sdev^2)) # vector of cumulative variance explained
  pcs <- which(cumvar >= threshold)[1]
  if(pcs<3) pcs <- 3

  ## standardize the pcs to the variance of the analog
  pcs.mean <- pcs.analog[, lapply(.SD, mean, na.rm = TRUE)]
  pcs.sd <- pcs.analog[, lapply(.SD, sd, na.rm = TRUE)]
  pcs.analog[, (names(pcs.analog)) := lapply(names(pcs.analog), function(col) {
    (get(col) - unlist(pcs.mean)[col]) / unlist(pcs.sd)[col]
  })]
  pcs.target[, (names(pcs.target)) := lapply(names(pcs.target), function(col) {
    (get(col) - unlist(pcs.mean)[col]) / unlist(pcs.sd)[col]
  })]

  ## Mahalanobis distance and sigma dissimilarity
  md <- (mahalanobis(pcs.target[,1:pcs], rep(0, pcs), var(pcs.analog[,1:pcs])))^0.5
  p <- pchi(md,pcs) # percentiles of the M distances on the chi distribution with degrees of freedom equaling the dimensionality of the distance measurement (PCs)
  q <- qchi(p,1) # values of the chi percentiles on a standard half-normal distribution (chi distribution with one degree of freedom)
  q[!is.finite(q)] <- 8 # set infinite values to 8 sigma (outside the decimal precision of pchi) 
  
  return(q)
  
  ## Plot
  if(plot2d){
    par(mfrow=c(2,2))
    for(i in 2:5){
    a <- pcs.analog[, c(1, i), with = FALSE]
    b <- pcs.target[, c(1, i), with = FALSE]
    plot(a, col="dodgerblue", xlim=range(c(a[,1], b[,1])), ylim=range(c(a[,2], b[,2])), asp=1)
    points(b, bg=ColScheme[cut(q, breakpoints)], pch=21, cex=1.5)
    points(a, col="dodgerblue", pch=16)
    mtext(paste(bgc.focal, "\n", pcs, "PCs"), line=-2.5, adj = 0.05, )
    }
  }
  
  ## 3D plot
  if(plot3d){
    
    # Extract the first three principal components
    a <- pcs.analog[, 1:3]
    b <- pcs.target[, 1:3]
    
    b_colors <- ColScheme[cut(q, breakpoints)] # Define colors for points in 'b'
    
    # Create the 3D scatterplot
    plot <- plot_ly() %>%
      add_trace(
        x = a[, 1], y = a[, 2], z = a[, 3],
        type = "scatter3d", mode = "markers",
        marker = list(size = 5, color = "dodgerblue", opacity = 1),
        name = "Analog Points"
      ) %>%
      add_trace(
        x = b[, 1], y = b[, 2], z = b[, 3],
        type = "scatter3d", mode = "markers",
        marker = list(size = 6, color = b_colors, opacity = 1),
        name = "Target Points"
      ) %>%
      layout(
        scene = list(
          xaxis = list(title = "PC1"),
          yaxis = list(title = "PC2"),
          zaxis = list(title = "PC3")
        ),
        title = list(text = bgc.focal, x = 0.05)
      )
    
    # Display the plot
    plot
  }
  
}

# test 
vars <- as.vector(outer(elements, c("wt", "sm"), paste, sep = "_"))
par(mar=c(1,1,1,1), mfrow=c(2,2))
for(i in c(2,3,4,6)){
  pcs <- i
  clim.analog <- clim.pts[pts$BGC==bgc.focal, ..vars]
  clim.target <- clim.grid[bgc.pred==bgc.focal, ..vars]
  sigma.target <- sigmaD(clim.target, clim.analog, pcs)
  pca <- prcomp(clim.analog, scale=TRUE)
  a <- predict(pca, clim.analog)[,1:2]
  b <- predict(pca, clim.target)[,1:2]
  plot(a, col="dodgerblue", xlim=range(c(a[,1], b[,1])), ylim=range(c(a[,2], b[,2])))
  points(b, bg=ColScheme[cut(sigma.target, breakpoints)], pch=21, cex=1.5)
  points(a, col="dodgerblue", pch=16)
  mtext(paste(bgc.focal, "\n", pcs, "PCs"), line=-2.5, adj = 0.05, )
}

par(mar=c(1,1,1,1), mfrow=c(1,2))
threshold = 0.95
for(bgc.focal in bgcs.target){
  clim.analog <- clim.pts[pts$BGC==bgc.focal, ..vars]
  clim.target <- clim.grid[bgc.pred==bgc.focal, ..vars]
  sigma[bgc.pred==bgc.focal] <- sigmaD_spattemp(clim.target, clim.analog, threshold)
}
X[clim.grid[, id]] <- sigma
plot(X, col=ColScheme, axes=F) 
mtext("Sigma novelty", side=4, line=-4.5, adj = 0.5, font=2)
mtext(paste0("spatiotemporal PCA method", "\n", "Threshold = ", threshold), line=-3.5, adj = 0.975)

pcs=4
for(bgc.focal in bgcs.target){
  clim.analog <- clim.pts[pts$BGC==bgc.focal, ..vars]
  clim.target <- clim.grid[bgc.pred==bgc.focal, ..vars]
  sigma[bgc.pred==bgc.focal] <- sigmaD(clim.target, clim.analog, pcs)
}
X[clim.grid[, id]] <- sigma
plot(X, col=ColScheme, axes=F) 
mtext("Sigma novelty", side=4, line=-4.5, adj = 0.5, font=2)
mtext(paste0("simple PCA method", "\n", "PCs = ", pcs), line=-3.5, adj = 0.975, )




#---------------------------
#---------------------------
# general function with multiple options
#---------------------------
#---------------------------

bgc.focal <- "CWHxm_WA" # good example of complete novelty
bgc.focal <- "CRFdh_CA" # good example of complete novelty
vars <- pred_vars

#' Novelty (outlyingness) measurement for pre-determined spatial climate analogs. 
#'
#' @description
#' 
#'
#' @details
#' 
#' This function optionally renders plots visualizing the climatic distance measurement in 
#' 2 and/or 3 dimensions.  
#' 
#' @param clim.targets data.table. Climate variables for which the analogs were 
#' identified
#' @param clim.analogs data.table. Climate variables of at least 50 locations 
#' representing the spatial variation in the reference period (historical) climate
#' of each analog in the analog pool. 
#' @param label.targets character. Vector of the analog IDs identified for the 
#' climatic conditions listed in `clim.targets`. Length equals number of records 
#' in `clim.targets`.  
#' @param label.analogs character. Vector of the analog IDs for the climatic 
#' conditions listed in `clim.analogs`. Length equals number of records in
#' `clim.analogs`.  
#' @param vars character. Climate variables to use in the novelty measurement. 
#' @param analog.focal character. Optionally specify a single analog for visualization. 
#' @param threshold numeric. The cumulative variance explained to use as a threshold
#' for truncation of principal components to use in the mahalanobis distance measurement. 
#' @param pcs integer. The fixed number of PCs to use in the mahalanobis distance 
#' measurement. Non-null values of this parameter override the `threshold` parameter
#' @param plot2d logical. If `analog.focal` is specified, plot a 4-panel set of bivariate plots visualizing the 
#' distances in the first five PCs. 
#' @param plot3d logical. If `analog.focal` is specified, plot a 3-dimensional plot 
#' visualizing the distances in the first three PCs. 
#'
#' @return `vector` of sigma dissimilarity (if sigma==TRUE) or Mahalanobis distances 
#' (if sigma==FALSE) corresponding to each element of the 'analogs.target' vector. 
#'
#' @examples
#' if (FALSE) {
#'   
#' }
#' #'
#' @export


clim.targets <- clim.grid
clim.analogs <- clim.pts
label.targets <- bgc.pred
label.analogs <- pts$BGC
vars <- pred_vars
analog.focal <- "IDFmw2"

analog_novelty <- function(clim.targets, clim.analogs, label.targets, label.analogs, vars, 
                           analog.focal = NULL, threshold = 0.95, pcs = NULL, 
                           plotScree = FALSE, plot2d = FALSE, plot3d = FALSE, 
                           plot2d.pcs = c(2,3,4,5), plot3d.pcs=c(1,2,3)){
  
  analogs <- if(is.null(analog.focal)) unique(label.targets) else analog.focal # list of analogs to loop through
  novelty <- rep(NA, length(label.targets)) # initiate a vector to store the sigma dissimilarities
  
  for(analog in analogs){ # loop through all of the analogs used to describe the target climates. 
    clim.analog <- clim.pts[label.analogs==analog, ..vars]
    clim.target <- clim.grid[label.targets==analog, ..vars]
    
    ## data cleaning
    clim.analog <- clim.analog[complete.cases(clim.analog)] # remove rows without data
    clim.analog <- clim.analog[, .SD, .SDcols = which(sapply(clim.analog, function(x) var(x, na.rm = TRUE) > 0))]  # Remove zero-variance columns
    clim.target <- clim.target[, .SD, .SDcols = names(clim.analog)]
    
    ## scale the data to the variance of the analog, since this is what we will ultimately be measuring the M distance in. 
    clim.mean <- clim.analog[, lapply(.SD, mean, na.rm = TRUE)]
    clim.sd <- clim.analog[, lapply(.SD, sd, na.rm = TRUE)]
    clim.analog[, (names(clim.analog)) := lapply(names(clim.analog), function(col) {
      (get(col) - unlist(clim.mean)[col]) / unlist(clim.sd)[col]
    })]
    clim.target[, (names(clim.target)) := lapply(names(clim.target), function(col) {
      (get(col) - unlist(clim.mean)[col]) / unlist(clim.sd)[col]
    })]
    
    ## PCA on pooled target and analog
    s <- sample(1:dim(clim.target)[1], dim(clim.analog)[1], replace = TRUE) # select a random sample of the target population to match the analog points. bootstrap if target population is smaller than analog points
    clim.target.sample <- clim.target[s,]
    pca <- prcomp(rbind(clim.analog, clim.target.sample), scale=FALSE)
    pcs.analog <- data.table(predict(pca, clim.analog))
    pcs.target <- data.table(predict(pca, clim.target))
    
    if(is.null(pcs)){
      ## select number of pcs
      cumvar <- cumsum(pca$sdev^2 / sum(pca$sdev^2)) # vector of cumulative variance explained
      pcs <- which(cumvar >= threshold)[1]
      if(pcs<3) pcs <- 3
    }
    
    ## z-standardize the pcs to the variance of the analog. this is necessary for a metric that can be translated into sigma values. 
    pcs.mean <- pcs.analog[, lapply(.SD, mean, na.rm = TRUE)]
    pcs.sd <- pcs.analog[, lapply(.SD, sd, na.rm = TRUE)]
    pcs.analog[, (names(pcs.analog)) := lapply(names(pcs.analog), function(col) {
      (get(col) - unlist(pcs.mean)[col]) / unlist(pcs.sd)[col]
    })]
    pcs.target[, (names(pcs.target)) := lapply(names(pcs.target), function(col) {
      (get(col) - unlist(pcs.mean)[col]) / unlist(pcs.sd)[col]
    })]
    
    ## Mahalanobis distance and sigma dissimilarity
    md <- (mahalanobis(pcs.target[,1:pcs], rep(0, pcs), var(pcs.analog[,1:pcs])))^0.5
    p <- pchi(md,pcs) # percentiles of the M distances on the chi distribution with degrees of freedom equaling the dimensionality of the distance measurement (PCs)
    q <- qchi(p,1) # values of the chi percentiles on a standard half-normal distribution (chi distribution with one degree of freedom)
    q[!is.finite(q)] <- 8 # set infinite values to 8 sigma (outside the decimal precision of pchi) 
    
    ## populate the novelty vector
    novelty[label.targets==analog] <- q
    
  } # end of the for-loop
  
  ## Plots for the final iteration of the for loop
  
  ## Scree plot
  if(plotScree){
    par(mfrow=c(1,1), mar=c(3,3,1,1), mgp=c(1.75,0.25,0))
    a <- apply(predict(pca, clim.analog), 2, sd)
    b <- apply(predict(pca, clim.target), 2, sd)
    diff <- abs(apply(predict(pca, clim.target), 2, mean) - apply(predict(pca, clim.analog), 2, mean))
    plot(0, xlim=c(1,length(a)), ylim=c(0,max(c(a,b, diff))*1.02), yaxs="i", col="white", tck=-0.005,
         xlab="Principal Component (PC)", ylab = "Standard Deviation")
    rect(pcs+0.5, -99, 99, 99, col = "grey95", lty=2)
    points(a, pch=21, bg="dodgerblue", cex=1.6)
    points(b, bg="grey", pch=21, cex=1.3)
    points(diff, col="black", pch=17, cex=1.3)
    text(pcs+0.5, max(c(a,b, diff)), paste0("Truncation at ", pcs, " PCs"), pos=4)
    legend("topright", legend=c("Analog", "Target", "Separation of means"), pt.bg=c("dodgerblue", "grey", NA), col = c("black", "black", "black"), pt.cex=c(1.6,1.3,1.3), pch=c(21, 21, 17), bty="n")
    box()
  }
  
  ## 2D scatterplot
  if(plot2d){
    par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(1.75,0.25,0))
    for(i in plot2d.pcs){
      a <- predict(pca, clim.analog)[, c(1, i)]
      b <- predict(pca, clim.target)[, c(1, i)]
      plot(a, col="dodgerblue", xlim=range(c(a[,1], b[,1])), ylim=range(c(a[,2], b[,2])), asp=1, tck=0.01)
      points(b, bg=ColScheme[cut(q, breakpoints)], pch=21, cex=1.5)
      points(a, col="dodgerblue", pch=16)
      mtext(paste(analog, "\n", pcs, "PCs"), line=-2.5, adj = 0.05, )
    }
  }
  
  ## 3D scatterplot
  if(plot3d){
    
    # revert to the raw pcs, because standardization obscures the shape of the analog distribution
    a <- predict(pca, clim.analog)
    b <- predict(pca, clim.target)
    
    b_colors <- ColScheme[cut(q, breakpoints)] # Define colors for points in 'b'
    
    # Create the 3D scatterplot
    plot <- plot_ly() %>%
      add_trace(
        x = a[, plot3d.pcs[1]], y = a[, plot3d.pcs[2]], z = a[, plot3d.pcs[3]],
        type = "scatter3d", mode = "markers",
        marker = list(size = 5, color = "dodgerblue", opacity = 1),
        name = "Analog Points"
      ) %>%
      add_trace(
        x = b[, plot3d.pcs[1]], y = b[, plot3d.pcs[2]], z = b[, plot3d.pcs[3]],
        type = "scatter3d", mode = "markers",
        marker = list(size = 6, color = b_colors, opacity = 1),
        name = "Target Points"
      ) %>%
      layout(
        scene = list(
          xaxis = list(title = paste0("PC", plot3d.pcs[1])),
          yaxis = list(title = paste0("PC", plot3d.pcs[2])),
          zaxis = list(title = paste0("PC", plot3d.pcs[3]))
        ),
        title = list(text = paste(analog,"\nNovelty in", pcs, "PCs"), x = 0.05)
      )
    
    # Display the plot
    print(plot)
  }

  return(novelty)
  
}

# plots of focal analogs
pcnum <- 4
bgc.focal = "CWHxm_WA" # moderate to high novelty
bgc.focal <- "CRFdh_CA"
bgc.focal <- "ESSFwm3" # good example of true analog plus novel fringe
bgc.focal <- "SBAPcp" 
bgc.focal <- "MGPmg" 
bgc.focal <- "BWBScmE" # good example of novelty along a low-variance dimension. 
bgc.focal <- "IDFmw2" # significant separation in the 4th PC. 
analog_novelty(clim.targets = clim.grid, 
               clim.analogs = clim.pts, 
               label.targets = bgc.pred, 
               label.analogs = pts$BGC, 
               vars = pred_vars, 
               pcs = pcnum,
               analog.focal = bgc.focal,
               plotScree = TRUE, 
               plot2d = TRUE,
               plot3d = TRUE,
               plot3d.pcs=c(1,2,4)
)


# Maps at increasing PCs
par(mar=c(1,1,1,1), mfrow=c(2,2))
pcnums=c(2,3,4,6)
for(pcnum in pcnums){
  novelty <- analog_novelty(clim.targets = clim.grid, 
                            clim.analogs = clim.pts, 
                            label.targets = bgc.pred, 
                            label.analogs = pts$BGC, 
                            vars = pred_vars, 
                            pcs = pcnum
                            
  )
  X[clim.grid[, id]] <- novelty
  plot(X, col=ColScheme, axes=F) 
  mtext("Sigma novelty", side=4, line=-4.5, adj = 0.5, font=2)
  mtext(paste0("RF variables", "\n", pcnum, "PCs"), line=-3.5, adj = 0.975, )
}

# Maps at defined PC
par(mar=c(1,1,1,1), mfrow=c(1,1))
pcnum=3
  novelty <- analog_novelty(clim.targets = clim.grid, 
                            clim.analogs = clim.pts, 
                            label.targets = bgc.pred, 
                            label.analogs = pts$BGC, 
                            vars = pred_vars, 
                            pcs = pcnum
                            
  )
  X[clim.grid[, id]] <- novelty
  plot(X, col=ColScheme, axes=F) 
  mtext("Sigma novelty", side=4, line=-4.5, adj = 0.5, font=2)
  mtext(paste0("RF variables", "\n", pcnum, "PCs"), line=-3.5, adj = 0.975, )

