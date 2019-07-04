################################################################################################################
# This R-script is used to compare VesselGraphs (see vesselgraph.h) datasets by visualizing the distribution
# of properties as density functions and boxplot diagrams.
#
# All input files are expected to be named "edges.csv" and "nodes.csv". The arguments supplied must be folders
# that contain such csv files. The basename of the folder is used to name the dataset.
################################################################################################################
separateFileExport = FALSE
# adapted from http://stackoverflow.com/questions/1366853/how-to-superimpose-multiple-density-curves-into-one-plot-in-r
compare <- function(data, stat_output_file, linenames, columnname, colors, main, xlab)
{

    s = list()
    alldata = NULL
    #Collect data
    for(i in 1:length(data))
    {
        d = data[[i]][[columnname]]
        s[[i]] = d
        alldata = c(d, alldata)
    }

    #Determine final bandwidth
    valuerange = range(alldata)
    bw = (valuerange[2] - valuerange[1]) / 32
    #bw = "nrd0"

    #Determine range for drawing
    junk.x = NULL
    junk.y = NULL
    for(i in 1:length(s))
    {
        junk.x = c(junk.x, density(s[[i]], bw = bw)$x)
        junk.y = c(junk.y, density(s[[i]], bw = bw)$y)
    }
    xr <- range(junk.x)
    yr <- range(junk.y)

    sep <- "================================================================================\n"
    cat(sep, file = stat_output_file, append = TRUE)
    cat(paste("Median(", main, "), ", xlab, "\n", sep=''), file = stat_output_file, append = TRUE)
    cat(sep, file = stat_output_file, append = TRUE)
    #Draw and print
    if(separateFileExport) {
        pdf(paste(columnname, "_density.pdf", sep=""))
    }
    plot(density(s[[1]], bw = bw), xlim = xr, ylim = yr, main = main, xlab = xlab)
    for(i in 1:length(s))
    {
        cat(paste(linenames[[i]], ": ", median(s[[i]]), "\n", sep=''), file = stat_output_file, append = TRUE)
        lines(density(s[[i]], bw = bw), xlim = xr, ylim = yr, col = colors[[i]])
        #abline(v = mean(s[[i]]), col = colors[[i]])
    }
    legend("topright", pch = "-", col = colors, legend = linenames)
    if(separateFileExport) {
        dev.off()
    }
    cat(paste(sep, "\n", sep=""), file = stat_output_file, append = TRUE)

    if(separateFileExport) {
        pdf(paste(columnname, "_boxplot.pdf", sep=""))
    }
    boxplot(s, col = colors, names = linenames, main = main, ylab = xlab)
    if(separateFileExport) {
        dev.off()
    }
}

plot.multi.bars.relative <- function(data, linenames, columnname, colors, main)
{

    s = list()
    alldata = NULL
    #Collect data
    max = 0
    for(i in 1:length(data)) {
        max = max(max, max(data[[i]][[columnname]]))
    }
    for(i in 1:length(data))
    {
        # Dirty hack: tabulate misses zeros, so we add 1 to all values
        d = tabulate(data[[i]][[columnname]]+1, max+1)
        d = d/NROW(data[[i]])
        s[[i]] = d
        alldata = c(d, alldata)
    }
    s = matrix(unlist(s), ncol=(max+1), byrow=TRUE)
    rownames(s) = linenames
    colnames(s) = 0:max
    #Draw
    if(separateFileExport) {
        pdf(paste(columnname, "_bars.pdf", sep=""))
    }
    barplot(s, main = main, col = colors, beside=TRUE, xlab = main, legend = linenames)
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

    # Remove one extreme outlier
    # n_data <- n_data[n_data$volume < 0.003,]
    return (n_data)
}
normalize_node_data <- function(data) {
    # Remove nodes at the border at the volume
    n_data <- data[data$isAtSampleBorder == 0,]
    return (n_data)
}
edge_data_file_name = "edges.csv"
node_data_file_name = "nodes.csv"

args <- commandArgs(trailingOnly = TRUE)
if(length(args)<1){
    stop("Expected parameters (i.e. following --args): [data_dir ]*")
}
colors <- 2:(length(args)+1)
data <- list()
node_data <- list()
names <- list()
for(i in 1:length(args))
{
    # Read and normalize edge data
    newData <- read.csv(file=file.path(args[[i]], edge_data_file_name), header=TRUE, sep=";")
    data[[i]] = normalize_edge_data(newData)
    print(args[[i]])
    print("edges")
    print(nrow(data[[i]]))
    # Calculate additional properties
    data[[i]]$straightness = 1/data[[i]]$curveness

    # Read and normalize node data
    newNodeData <- read.csv(file=file.path(args[[i]], node_data_file_name), header=TRUE, sep=";")
    node_data[[i]] = normalize_node_data(newNodeData)
    print("nodes")
    print(nrow(node_data[[i]]))

    names[[i]] = basename(basename(args[[i]])) #This is really not universal...
}

# Prepare output stat file
stat_output_file <- "stats.txt"
if (file.exists(stat_output_file)) file.remove(stat_output_file)

compare(data, stat_output_file = stat_output_file, linenames = names, columnname = "length", col = colors, main = "Length", xlab = "l/mm")

compare(data, stat_output_file = stat_output_file, linenames = names, columnname = "distance", col = colors, main = "Distance", xlab = "l/mm")

compare(data, stat_output_file = stat_output_file, linenames = names, columnname = "curveness", col = colors, main = "Curveness", xlab = "Curveness")
compare(data, stat_output_file = stat_output_file, linenames = names, columnname = "straightness", col = colors, main = "Straightness", xlab = "Straightness")


## Vessel crosssection related
# Avg
compare(data, stat_output_file = stat_output_file, linenames = names, columnname = "minRadiusAvg", col = colors, main = "Min Radius Mean", xlab = "l/mm")

compare(data, stat_output_file = stat_output_file, linenames = names, columnname = "avgRadiusAvg", col = colors, main = "Average Radius Mean", xlab = "l/mm")

compare(data, stat_output_file = stat_output_file, linenames = names, columnname = "maxRadiusAvg", col = colors, main = "Max Radius Mean", xlab = "l/mm")

#TODO 6?
#dat <- list(edge_dataa$minRadiusAvg, edge_data$avgRadiusAvg, edge_dataa$maxRadiusAvg)
#compare(data, linenames = c("min", "avg", "max"), col = c("red", "orange", "yellow"), main = "Average Radii", xlab = "l/mm")

# STD
compare(data, stat_output_file = stat_output_file, linenames = names, columnname = "minRadiusStd", col = colors, main = "Min Radius Standard Deviation", xlab = "l/mm")

compare(data, stat_output_file = stat_output_file, linenames = names, columnname = "avgRadiusStd", col = colors, main = "Average Radius Standard Deviation", xlab = "l/mm")

compare(data, stat_output_file = stat_output_file, linenames = names, columnname = "maxRadiusStd", col = colors, main = "Max Radius Standard Deviation", xlab = "l/mm")

# Roundness Avg
compare(data, stat_output_file = stat_output_file, linenames = names, columnname = "roundnessAvg", col = colors, main = "Roundness Mean", xlab = "Roundness")

# Roundness Avg
compare(data, stat_output_file = stat_output_file, linenames = names, columnname = "roundnessStd", col = colors, main = "Roundness Standard Deviation", xlab = "Roundness")

# Volume
compare(data, stat_output_file = stat_output_file, linenames = names, columnname = "volume", col = colors, main = "Volume", xlab = expression(paste("l"^"3"*"/mm"^"3")))

# Avg Cross section
compare(data, stat_output_file = stat_output_file, linenames = names, columnname = "avgCrossSection", col = colors, main = "Average Cross Section", expression(paste("l"^"2"*"/mm"^"2")))

# Degree of nodes
plot.multi.bars.relative(data, linenames = names, columnname = "node1_degree", col = colors, main = "Lower per-edge Node Degree")
plot.multi.bars.relative(data, linenames = names, columnname = "node2_degree", col = colors, main = "Higher per-edge Node Degree")

# Degree of nodes
plot.multi.bars.relative(node_data, linenames = names, columnname = "degree", col = colors, main = "Node Degree")

# Close output stat file
#close(stat_output_file)
