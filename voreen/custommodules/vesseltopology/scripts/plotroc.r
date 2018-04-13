# adapted from http://stackoverflow.com/questions/1366853/how-to-superimpose-multiple-density-curves-into-one-plot-in-r
plot.multi.roc <- function(data, linenames, colors, main)
{
    sensitivity = list()
    one_minus_specificity = list()
    thresholds = list()
    jaccard = list()
    f1score = list()
    #Collect data
    for(i in 1:length(data))
    {
        one_minus_specificity[[i]] = data[[i]]$falsePositive / data[[i]]$numVoxels
        precision = data[[i]]$truePositive / (data[[i]]$truePositive+data[[i]]$falsePositive)
        sensitivity[[i]] = data[[i]][["sensitivity"]]
        thresholds[[i]] = data[[i]][["threshold"]]
        jaccard[[i]] = data[[i]][["jaccardIndex"]]
        #f1score[[i]] = 2*data[[i]][["sensitivity"]]*precision/(data[[i]][["sensitivity"]]+precision)
        #f1score[[i]] = data[[i]]$truePositive/(data[[i]]$truePositive+data[[i]]$falseNegative+data[[i]]$falsePositive)
    }

    #Draw
    #plot(one_minus_specificity[[1]], sensitivity[[1]], main = main, xlab="False Positive Rate", ylab="Sensitivity", lty=2, lwd=0, type="n", xlim=c(0,0.2), ylim=c(0.8,1))
    ##plot(one_minus_specificity[[1]], sensitivity[[1]], main = main, xlab="False Positive Rate", ylab="Sensitivity", lty=2, lwd=0, type="n")
    #for(i in 1:length(data))
    #{
    #    lines(one_minus_specificity[[i]], sensitivity[[i]], lty=1, lwd=1, col = colors[[i]])
    #    points(one_minus_specificity[[i]], sensitivity[[i]], cex=1, pch=20, col = colors[[i]])
    #    text(one_minus_specificity[[i]], sensitivity[[i]], col = colors[[i]], labels=thresholds[[i]], pos=4)
    #}
    #legend("bottomright", pch = "-", col = colors, legend = linenames)

    plot(jaccard[[1]], main = "", xlab="Threshold", xaxt="n", ylab="Jaccard Index", lty=2, lwd=0, type="n")
    for(i in 1:length(data))
    {
        lines(jaccard[[i]], lty=1, lwd=1, col = colors[[i]])
        axis(1, at=c(0,50), labels=c("low", "high"))
        #text(jaccard[[i]], col = colors[[i]], labels=thresholds[[i]], pos=4)
        print(paste(linenames[[i]], max(jaccard[[i]])))
    }
    legend("topright", pch = "-", col = colors, legend = linenames)
}

args <- commandArgs(trailingOnly = TRUE)
if(length(args)<1){
    stop("Expected parameters (i.e. following --args): [data_dir ]*")
}
data <- list()
node_data <- list()
names <- list()
for(i in 1:length(args))
{
    newData <- read.csv(file=args[[i]], header=TRUE, sep=";")
    data[[i]] = newData
    names[[i]] = (sub("^([^.]*).*", "\\1", basename(args[[i]])))
}
#colors <- c("red", "green", "blue", "orange", "black", "grey")
colors <- 1:length(args)

plot.multi.roc(data, linenames = names, col = colors, main = "ROC")
