################################################################################################################
# This R-script performs a wilcox test to decide whether there is a statistical significant difference
# between distribution of edge properties in vessel graphs.
#
# Paths to the edges.csv datasets that are to be compared are expected as parameters.
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
if(length(args)!=2){
    stop("Expected parameters: dataset1 dataset2")
}
alldata1 = normalize_edge_data(read.csv(file=file.path(args[[1]]), header=TRUE, sep=";"))
alldata2 = normalize_edge_data(read.csv(file=file.path(args[[2]]), header=TRUE, sep=";"))

alldata1$straightness = 1/alldata1$curveness
alldata2$straightness = 1/alldata2$curveness


props = c('length','distance','curveness','straightness','volume','avgCrossSection','minRadiusAvg','minRadiusStd','avgRadiusAvg','avgRadiusStd','maxRadiusAvg','maxRadiusStd','roundnessAvg','roundnessStd')

print("Performing Mann–Whitney U test:")

for (prop in props) {
    dataset1 = alldata1[[prop]]
    dataset2 = alldata2[[prop]]

    test = wilcox.test(dataset1, dataset2, exact=TRUE)
    p_value = test$p.value*length(props) #Bonferroni correction
    s = sprintf('%s: p=%e', prop, p_value)
    print(s)
}
