#include "lib.h"
#include "network.h"
#include <cuda.h>
#include "../../datastructures/vesselgraph.h"

#include <iostream>

stim::network<float> GT;			// ground truth network
stim::network<float> T;			// test network

int device = -1;

void setdevice(){
	int count;
	cudaGetDeviceCount(&count);					// numbers of device that are available
	if(count > 0) {
		device = 1;
	} else {
		std::cout<<"No cuda device available."<<std::endl;
		device = -1;
	}
}

float netmets_compare_networks(const voreen::VesselGraph& g1, const voreen::VesselGraph& g2) {
    if(device == -1) {
	setdevice();
    }
    // find appropriate radius
    float radiusSum = 0.0;
    for(const auto& edge : g1.getEdges()) {
        radiusSum += edge.getAvgRadiusAvg();
    }
    for(const auto& edge : g2.getEdges()) {
        radiusSum += edge.getAvgRadiusAvg();
    }
    float globalAvgRadius = radiusSum / (g1.getEdges().size() + g2.getEdges().size());
    float sigma = globalAvgRadius;

    GT = GT.compare(T, sigma, device);				// compare the ground truth to the test case - store errors in GT
T = T.compare(GT, sigma, device);				// compare the test case to the ground truth - store errors in T

    //calculate the metrics
    float FPR = GT.average();						// calculate the metrics
    float FNR = T.average();

    return FPR*FNR;
}
