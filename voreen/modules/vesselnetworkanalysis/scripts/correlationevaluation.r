################################################################################################################
# This R-script is used to evaluate/find/visualize correlations between properties of VesselGraphs
# (see vesselgraph.h) by creating scatter plots and a correlation matrix.
#
# All input files are expected to be named "edges.csv". The arguments supplied must be folders that contain
# such a csv file. The basename of the folder is used to name the dataset.
################################################################################################################
separateFileExport = FALSE

plot.multi.scatter <- function(data, datasetnames, xcolumnname, ycolumnname, prettyxname=xcolumnname, prettyyname=ycolumnname)
{
    s = list()
    #Collect data

    xdata <- NULL
    ydata <- NULL
    colordata <- NULL
    for(i in 1:length(data))
    {
        d = data[[i]]
        #print(NROW(d))

        xdata <- c(xdata, d[[xcolumnname]])
        ydata <- c(ydata, d[[ycolumnname]])
        colordata <- c(colordata, rep(i, each=nrow(ydata)))
    }

    if(separateFileExport) {
        pdf(paste(xcolumnname, "_vs_", ycolumnname,"_correlation.pdf", sep=""))
    }
    plot(xdata, ydata, col=colordata, xlab=prettyxname, ylab=prettyyname, main= " ")
    legend("topright", pch = 1, col = 1:length(data), legend = datasetnames)
    #abline(lm(ydata ~ xdata), col="red")
    #lines(lowess(xdata, ydata), col="blue") # lowess line (x,y)
    if(separateFileExport) {
        dev.off()
    }
}

plot.correlations <- function(data, columnnames, prettycolumnnames = columnnames)
{
    n = length(columnnames)
    corrMat = matrix( nrow=n, ncol=n)
    method = "pearson"
    for(i in 1:n)
    {
        for(j in 1:n)
        {
            x <- NULL
            y <- NULL
            for(k in 1:length(data))
            {
                x = c(x, data[[k]][[columnnames[[i]]]])
                y = c(y, data[[k]][[columnnames[[j]]]])
            }
            corrMat[i,j] = cor(x,y, method=method)
        }
    }
    if(separateFileExport) {
        pdf("correlation_matrix.pdf")
    }

    par(mar = c(7, 10, 1, 1))
    colorPalette = colorRampPalette(c("red", "white", "blue"))(101)
    image(t(corrMat),col=colorPalette, zlim=c(-1,1), x=1:n, y=1:n, xlab=" ", ylab=" ", xaxt="n", yaxt="n")
    print(corrMat)

    for(i in 1:n)
    {
        for(j in 1:n)
        {
            x <- NULL
            y <- NULL
            for(k in 1:length(data))
            {
                x = c(x, data[[k]][[columnnames[[i]]]])
                y = c(y, data[[k]][[columnnames[[j]]]])
            }
            val = cor(x,y, method = method)
            text(i,j, sprintf("%0.2f", val), col = "black")
        }
    }

    axis(1, at=1:n, labels = FALSE)
    text(1:n, par("usr")[1] - 0.25, srt = 45, adj = 1,
         labels = prettycolumnnames, xpd = TRUE, cex=0.65)

    axis(2, at=1:n, labels = FALSE)
    text(par("usr")[1] - 0.25, 1:n, srt = 0, adj = 1,
         labels = prettycolumnnames, xpd = TRUE, cex=0.65)
    if(separateFileExport) {
        dev.off()
    }
}
normalize_edge_data <- function(data) {
    n_data <- data
    # Remove very short edges
    n_data <- n_data[n_data$num_voxels > 5,]
    # Remove edges corresponding to small "blobs"
    n_data <- n_data[n_data$length > 2*n_data$avgRadiusAvg,]
    # Remove loops
    n_data <- n_data[n_data$distance != 0,]
    # Remove edges that do not have radius data
    n_data <- n_data[n_data$maxRadiusAvg != -1,]
    # Remove edges that leave the volume
    n_data <- n_data[n_data$hasNodeAtSampleBorder == 0,]
    return (n_data)
}

edge_data_file_name = "edges.csv"

args <- commandArgs(trailingOnly = TRUE)
if(length(args)<1){
    stop("Expected parameters (i.e. following --args): [data_dir ]*")
}
data <- list()
names <- list()
for(i in 1:length(args))
{
    newData <- read.csv(file=file.path(args[[i]], edge_data_file_name), header=TRUE, sep=";")
    newData <- normalize_edge_data(newData)
    data[[i]] = newData
    names[[i]] = basename(args[[i]])

    # Calculate additional properties
    data[[i]]$straightness = 1/data[[i]]$curveness
}

#length group
plot.multi.scatter(data, names, "length", "distance")
plot.multi.scatter(data, names, "length", "curveness")
plot.multi.scatter(data, names, "distance", "curveness")
plot.multi.scatter(data, names, "length", "volume")

plot.multi.scatter(data, names, "length", "roundnessAvg", "length", "meanroundness")
plot.multi.scatter(data, names, "length", "avgRadiusStd", "length", "stdroundness")

#roundness group
plot.multi.scatter(data, names, "roundnessAvg", "roundnessStd", "meanroundness", "stdroundness")
plot.multi.scatter(data, names, "roundnessAvg", "avgRadiusAvg", "meanroundness", expression(meanrad[avg]))
plot.multi.scatter(data, names, "roundnessAvg", "maxRadiusAvg", "meanroundness", expression(meanrad[max]))
plot.multi.scatter(data, names, "avgRadiusAvg", "maxRadiusAvg", expression(meanrad[avg]), expression(meanrad[max]))
plot.multi.scatter(data, names, "avgRadiusAvg", "avgCrossSection", expression(meanrad[avg]), "cross-section")
plot.multi.scatter(data, names, "roundnessAvg", "avgCrossSection", "meanroundness", "cross-section")

plot.multi.scatter(data, names, "roundnessAvg", "avgRadiusStd", "meanroundness", expression(stdrad[avg]))
plot.multi.scatter(data, names, "roundnessAvg", "maxRadiusStd", "meanroundness", expression(stdrad[max]))
plot.multi.scatter(data, names, "avgRadiusAvg", "avgRadiusStd", expression(meanrad[avg]), expression(stdrad[avg]))
plot.multi.scatter(data, names, "maxRadiusAvg", "maxRadiusStd", expression(meanrad[max]), expression(stdrad[max]))

#radius std group
plot.multi.scatter(data, names, "avgRadiusStd", "maxRadiusStd", expression(stdrad[avg]), expression(stdrad[max]))

rowNames = c( "length", "distance", "straightness", "roundnessAvg", "roundnessStd", "volume", "avgCrossSection", "avgRadiusAvg", "maxRadiusAvg", "avgRadiusStd", "maxRadiusStd")
#prettynames = c(
#                "length",
#                "distance",
#                "straightness",
#                "meanroundness",
#                "stdroundness",
#                "volume",
#                "cross-section",
#                expression(meanrad[avg]),
#                expression(meanrad[max]),
#                expression(stdrad[avg]),
#                expression(stdrad[max])
#                          )
prettynames = c(
                "segment length",
                "distance",
                "straightness",
                "roundness mean",
                "roundness standard deviation",
                "volume per segment",
                "average cross-section",
                "average radius mean",
                "maximal radius mean",
                "average radius standard deviation",
                "maximal radius standard deviation"
                          )
plot.correlations(data, rowNames, prettynames)
