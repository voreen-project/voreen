################################################################################################################
# This R-script is used to visualize properties of VesselGraphs (see vesselgraph.h).
#
# The first argument must be the path to the edge data. The second must be the path to the node data.
################################################################################################################
# adapted from http://stackoverflow.com/questions/1366853/how-to-superimpose-multiple-density-curves-into-one-plot-in-r
plot.multi.dens <- function(s, linenames, col, main, xlab)
{
    junk.x = NULL
    junk.y = NULL
    for(i in 1:length(s))
    {
        junk.x = c(junk.x, density(s[[i]])$x)
        junk.y = c(junk.y, density(s[[i]])$y)
    }
    xr <- range(junk.x)
    yr <- range(junk.y)
    plot(density(s[[1]]), xlim = xr, ylim = yr, main = main, xlab = xlab)
    for(i in 1:length(s))
    {
        lines(density(s[[i]]), xlim = xr, ylim = yr, col = col[[i]])
    }
    legend("topright", pch = "-", col = col, legend = linenames)
}

args <- commandArgs(trailingOnly = TRUE)
if(length(args)<2){
    stop("Expected parameters (i.e. following --args): edge_data.csv node_data.csv")
}

edge_file = args[[1]]
node_file = args[[2]]

edge_data <- read.csv(file=edge_file, header=TRUE, sep=";")

# Remove very short edges
edge_data <- edge_data[edge_data$num_voxels > 2,]
# Remove loops
edge_data <- edge_data[edge_data$distance != 0,]


## Vessel length related
names(edge_data)
hist(edge_data$length, breaks=20, freq=FALSE, xlab="l/mm(?)", main="Average Length")
lines(density(edge_data$length), col="red")

hist(edge_data$distance, breaks=20, freq=FALSE, xlab="l/mm(?)", main="Average Distance")
lines(density(edge_data$distance), col="red")

plot(density(edge_data$distance), col="red", xlab="l/mm(?)", main="Average Distance vs Average Length")
lines(density(edge_data$length), col="green")
legend("topright", pch = "-", col = c("red", "green"), legend = c("Distance", "Length"))

hist(edge_data$curveness, breaks=20, freq=FALSE, xlab="Curveness", main="Average Curveness (0 => Loop!)")
lines(density(edge_data$curveness), col="red")


## Vessel crosssection related
# Avg
hist(edge_data$minRadiusAvg, freq=FALSE, xlab="l/mm", main="Average Min Radius")
lines(density(edge_data$minRadiusAvg), col="red")
hist(edge_data$avgRadiusAvg, freq=FALSE, xlab="l/mm", main="Average Average Radius")
lines(density(edge_data$avgRadiusAvg), col="red")
hist(edge_data$maxRadiusAvg, freq=FALSE, xlab="l/mm", main="Average Max Radius")
lines(density(edge_data$maxRadiusAvg), col="red")

dat <- list(edge_data$minRadiusAvg, edge_data$avgRadiusAvg, edge_data$maxRadiusAvg)
plot.multi.dens(dat, linenames = c("min", "avg", "max"), col = c("red", "orange", "yellow"), main = "Average Radii", xlab = "l/mm")

# STD
hist(edge_data$minRadiusStd, freq=FALSE, xlab="l/mm", main="Min Radius Standard Deviation")
lines(density(edge_data$minRadiusStd), col="red")
hist(edge_data$avgRadiusStd, freq=FALSE, xlab="l/mm", main="Average Radius Standard Deviation")
lines(density(edge_data$avgRadiusStd), col="red")
hist(edge_data$maxRadiusStd, freq=FALSE, xlab="l/mm", main="Max Radius Standard Deviation")
lines(density(edge_data$maxRadiusStd), col="red")

dat <- list(edge_data$minRadiusStd, edge_data$avgRadiusStd, edge_data$maxRadiusStd)
plot.multi.dens(dat, linenames = c("min", "avg", "max"), col = c("red", "orange", "yellow"), main = "Radii Standard Devation", xlab = "l/mm")

# Roundness Avg
hist(edge_data$roundnessAvg, freq=FALSE, xlab="Roundness", main="Average Roundness")
lines(density(edge_data$roundnessAvg), col="red")

# Roundness Std
hist(edge_data$roundnessStd, freq=FALSE, xlab="Roundness", main="Roundness Standard Deviation")
lines(density(edge_data$roundnessStd), col="red")

# comparison: roundness -> avg vs avg -> roundness
roundnessAvg <- edge_data$minRadiusAvg/edge_data$maxRadiusAvg
dat <- list(edge_data$roundnessAvg, roundnessAvg)
plot.multi.dens(dat, linenames = c("roundness -> avg", "avg -> roundness"), col = c("green", "red"), main = "roundness -> avg vs. avg -> roundness", xlab = "Roundness")

## Various Scatter plots
# Roundness Std<->Avg
xdata <- edge_data$roundnessAvg
xlabel <- "Average Roundness"
ydata <- edge_data$roundnessStd
ylabel <- "Roundness Standard Deviation"
plot(xdata, ydata, xlab=xlabel, ylab=ylabel, main="Roundness Scattered")
abline(lm(ydata ~ xdata), col="red")
lines(lowess(xdata, ydata), col="blue") # lowess line (x,y)

xdata <- edge_data$roundnessAvg
xlabel <- "Average Roundness"
ydata <- edge_data$avgRadiusStd
ylabel <- "Average Radius Standard Deviation"
plot(xdata, ydata, xlab=xlabel, ylab=ylabel, main="Roundness Scattered")
abline(lm(ydata ~ xdata), col="red")
lines(lowess(xdata, ydata), col="blue") # lowess line (x,y)

xdata <- edge_data$roundnessAvg
xlabel <- "Average Roundness"
ydata <- edge_data$curveness
ylabel <- "Curveness"
plot(xdata, ydata, xlab=xlabel, ylab=ylabel, main="Roundness Scattered")
abline(lm(ydata ~ xdata), col="red")
lines(lowess(xdata, ydata), col="blue") # lowess line (x,y)

xdata <- edge_data$avgRadiusStd
xlabel <- "Average Radius Standard Deviation"
ydata <- edge_data$curveness
ylabel <- "Curveness"
plot(xdata, ydata, xlab=xlabel, ylab=ylabel, main="Roundness Scattered")
abline(lm(ydata ~ xdata), col="red")
lines(lowess(xdata, ydata), col="blue") # lowess line (x,y)

# nodes
node_data <- read.csv(file=node_file, header=TRUE, sep=";")
degree_hist = table(node_data$degree)
plt = barplot(degree_hist, main = "Node Degree")
text(plt, degree_hist-1, paste("n = ",degree_hist,sep="") ,cex=1)
