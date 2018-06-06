################################################################################################################
# This R-script is used to compare VesselGraphs (see vesselgraph.h) datasets by visualizing the distribution
# of properties as density functions and boxplot diagrams.
#
# All input files are expected to be named "edges.csv" and "nodes.csv". The arguments supplied must be folders
# that contain such csv files. The basename of the folder is used to name the dataset.
################################################################################################################
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

    # Remove one extreme outlier
    # n_data <- n_data[n_data$volume < 0.003,]
    return (n_data)
}

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=3){
    stop("Expected parameters: dataset1 dataset2 property")
}
property = args[[3]]
alldata1 = normalize_edge_data(read.csv(file=file.path(args[[1]]), header=TRUE, sep=";"))
alldata2 = normalize_edge_data(read.csv(file=file.path(args[[2]]), header=TRUE, sep=";"))

alldata1$straightness = 1/alldata1$curveness
alldata2$straightness = 1/alldata2$curveness


dataset1 = alldata1[[property]]
dataset2 = alldata2[[property]]


shapiro.test(dataset1)
shapiro.test(dataset2)
#t.test(dataset1, dataset2, var.equal=FALSE, conf.level=0.95)

wilcox.test(dataset1, dataset2)
