################################################################################################################
# This R-script is used to visualize and compare the performance of different segmentations. The input file
# can be produced by a SegmentationListEvaluation processor (see segmentationlistevaluation.h)
#
# Path to the dataset files (in csv format) must be supplied as arguments
################################################################################################################
plot.multi.hist <- function(data, linenames, columnname, colors, main, xlab)
{
    s = list()
    #Collect data
    for(i in 1:length(data))
    {
        d = data[[i]][[columnname]]
        s[[i]] = d/255*0.10048846 #Hardcoded hack for specific dataset!
        #s[[i]] = d
    }

    #Draw
    for(i in 1:length(s))
    {
        #h = hist(s[[i]], main = main, xlab = xlab, ylab=" ")
        h = hist(s[[i]], plot=FALSE, breaks=64)
        print(h)
        barplot(pmax(log(h$counts),0), names.arg = h$mids, xlab = "Vesselness", ylab = "Log10(frequency)")
        #hist(s[[i]])
        #abline(v = mean(s[[i]]), col = colors[[i]])
    }
}

args <- commandArgs(trailingOnly = TRUE)
if(length(args)<1){
    stop("Expected parameters (i.e. following --args): [data_dir ]*")
}
colors <- c("red", "green", "blue")
data <- list()
node_data <- list()
names <- list()
for(i in 1:length(args))
{
    newData <- read.csv(file=args[[i]], header=TRUE, sep=";")
    data[[i]] = newData
    names[[i]] = basename(basename(args[[i]])) #This is really not universal...
}

plot.multi.hist(data, linenames = names, columnname = "value", col = colors, main = "Intensity", xlab = "Intensity")
