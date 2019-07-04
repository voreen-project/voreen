#include "lib.h"
#include "network.h"
#include <cuda.h>
#include "../../datastructures/vesselgraph.h"

#include <iostream>

int netmetsCudaDevice = -1;

void setdevice(){
    int count;
    cudaGetDeviceCount(&count);                    // numbers of device that are available
    if(count > 0) {
        netmetsCudaDevice = 1;
    } else {
        std::cout<<"No cuda device available."<<std::endl;
        netmetsCudaDevice = -1;
    }
}

typedef stim::network<float> Network;

Network networkFromVesselGraph(const voreen::VesselGraph& g) {
    Network output;

    for(const auto& node : g.getNodes()) {
        if(node.getDegree() == 0) {
            // Freestanding nodes are not handled well by netmets
            continue;
        }

        stim::vec3<float> pos(node.pos_.x, node.pos_.y, node.pos_.z);
        Network::vertex new_vertex(pos);

        for(const auto& edge : g.getEdges()) {
            if(edge.getNodeID1() == node.getID()) {
                new_vertex.e[0].push_back(edge.getID().raw());
            }
            if(edge.getNodeID2() == node.getID()) {
                new_vertex.e[1].push_back(edge.getID().raw());
            }
        }
        output.V.push_back(new_vertex);
    }

    for(const auto& edge : g.getEdges()) {
        const auto& voxels = edge.getVoxels();
        size_t centerline_length = voxels.size() + 2;
        stim::centerline<float> c3(centerline_length);
        std::vector<float> radius(c3.size());

        c3[0] = stim::vec3<float>(edge.getNode1().pos_.x, edge.getNode1().pos_.y, edge.getNode1().pos_.z);
        radius[0] = edge.getNode1().getRadius(); //TODO
        radius[0] = 0.0;
        int i = 1;
        for(const auto& voxel : voxels) {
            stim::vec3<float> pos(voxel.pos_.x, voxel.pos_.y, voxel.pos_.z);
            c3[i] = pos;
            if(voxel.hasValidData()) {
                radius[i] = voxel.avgDistToSurface_;
            } else {
                radius[i] = 0.0; //TODO maybe think of something else
            }
            ++i;
        }
        c3[centerline_length - 1] = stim::vec3<float>(edge.getNode2().pos_.x, edge.getNode2().pos_.y, edge.getNode2().pos_.z);
        radius[centerline_length - 1] = edge.getNode2().getRadius();
        radius[centerline_length - 1] = 0.0;

        c3.update();

        stim::cylinder<float> C3(c3);
        C3.copy_r(radius);

        Network::edge new_edge(C3);
        new_edge.v[0] = edge.getNodeID1().raw();
        new_edge.v[1] = edge.getNodeID2().raw();

        output.E.push_back(new_edge);
    }

    return output;
}

NetmetsResult netmets_compare_networks(const voreen::VesselGraph& groundtruthNetwork, const voreen::VesselGraph& testNetwork) {
    if(netmetsCudaDevice == -1) {
        setdevice();
    }
    // find appropriate radius
    float radiusSum = 0.0;
    size_t numConsideredVoxels = 0;
    auto processGraph = [&] (const voreen::VesselGraph& g) {
        for(const auto& edge : g.getEdges()) {
            for(const auto& voxel : edge.getVoxels()) {
                if(voxel.hasValidData()) {
                    radiusSum += voxel.avgDistToSurface_;
                    ++numConsideredVoxels;
                }
            }
        }
    };
    processGraph(groundtruthNetwork);
    float globalAvgRadius = radiusSum / numConsideredVoxels;
    float sigma = globalAvgRadius * 10.0 /*Found experimentally to yield good results in netmets application */;
    std::cout << "Using sigma: " << sigma << std::endl;

    auto gtnm = networkFromVesselGraph(groundtruthNetwork);
    auto tnm = networkFromVesselGraph(testNetwork);

    // For debugging purposes:
    //gtnm.saveNwt("/home/dominik/g1.nwt");
    //tnm.saveNwt("/home/dominik/g2.nwt");

    // TODO check that overwriting is fine
    gtnm = gtnm.compare(tnm, sigma, netmetsCudaDevice);                // compare the ground truth to the test case - store errors in GT
    tnm = tnm.compare(gtnm, sigma, netmetsCudaDevice);                // compare the test case to the ground truth - store errors in T

    NetmetsResult result;
    //calculate the metrics TODO: Make sure that fnr/fpr are in the right order here
    result.fnr = gtnm.average();
    result.fpr = tnm.average();

    return result;
}
