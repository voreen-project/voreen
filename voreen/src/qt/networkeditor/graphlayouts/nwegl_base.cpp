/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#include "voreen/qt/networkeditor/graphlayouts/nwegl_base.h"

#include "voreen/core/network/networkgraph.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/ports/port.h"
#include "voreen/core/processors/processor.h"

#include <set>
#include <queue>
#include <map>
#include <iostream>
#include <algorithm>

namespace voreen{

std::vector<Processor*> NWEGL_Base::getSourceProcessors(ProcessorNetwork* network) const {
    std::vector<Processor*> sourceProcessors;
    for (size_t i = 0; i <network->getProcessors().size(); ++i) {
        if (network->getProcessors()[i]->isSource() == true)
            sourceProcessors.push_back(network->getProcessors()[i]);
    }
    return sourceProcessors;
}

std::vector<Processor*> NWEGL_Base::getEndProcessors(ProcessorNetwork* network) const {
    std::vector<Processor*> endProcessors;
    for (size_t i = 0; i < network->getProcessors().size(); ++i) {
        if (network->getProcessors()[i]->isEndProcessor() == true)
            endProcessors.push_back(network->getProcessors()[i]);
    }
    return endProcessors;
}

std::vector<std::vector<Processor*> > NWEGL_Base::getGraphLayer(std::vector<Processor*> proc) const {

    //vector that saves the layers ordered from bottom to top
    std::vector<std::vector<Processor*> > layerList;

    //helping variables to determine predecessor processors
    std::vector<Port*> ports;
    std::vector<const Port*> connPort;
    std::vector<Processor*> connProc;

    //endProcesors for initialisation of the sorting algorithm
    std::vector<Processor*> endProc; /* = getEndProcessors(network);*/

    //find endProcessors
    for (size_t i = 0; i < proc.size(); ++i) {
        if (proc[i]->isEndProcessor() == true) {
            endProc.push_back(proc[i]);
        }
    }

    //if there are no endprocessors in the network (or none selected) take those which connected processors are not selected
    if (endProc.empty()) {
        for (size_t i = 0; i < proc.size(); ++i) {
            std::vector<Processor*> connProc;

            std::vector<Port*> out = proc[i]->getOutports();
            for (size_t k = 0; k < out.size(); k++) {
                // Ignore loop ports.
                if(!out[k]->isLoopPort()) {
                    std::vector<const Port *> in = out[k]->getConnected();
                    for (size_t l = 0; l < in.size(); l++) {
                        Processor *connect = in[l]->getProcessor();
                        connProc.push_back(connect);
                    }
                }
            }

            bool isSelected = false;
            for (size_t j = 0; j < connProc.size(); j++) {
                if(std::find(proc.begin(), proc.end(), connProc[j])!=proc.end()) {
                    isSelected = true;
                }
            }
            //make sure that coprocessors are not counted as endprocessors
            if (isSelected == false && proc[i]->getCoProcessorOutports().empty()) {
                endProc.push_back(proc[i]);
            }
        }
    }

    //at first save endProcessors in layerList
    layerList.push_back(endProc);

    //bool variable checks if all processors have a layer
    bool lastConn = false;

    //now add next layers to layerList
    while(lastConn == false) {

        //saves actual layer
        std::vector<Processor*> actLayer = layerList.back();

        for (size_t layerIt = 0; layerIt < actLayer.size(); layerIt++) {
            ports = actLayer[layerIt]->getInports();
            for (size_t portIt = 0; portIt < ports.size(); portIt++) {
                if(!ports[portIt]->isLoopPort()) {
                    connPort = ports[portIt]->getConnected();
                    for (size_t procIt = 0; procIt < connPort.size(); procIt++) {
                        //erase if same processor is inside layer multiple times (and only if processor is a selected one)
                        if (std::find(proc.begin(), proc.end(), connPort[procIt]->getProcessor()) != proc.end()) {
                            for (size_t it = 0; it < connProc.size(); it++) {
                                if (connPort[procIt]->getProcessor() == connProc[it]) {
                                    connProc.erase(connProc.begin() + it);
                                }
                            }
                            connProc.push_back(connPort[procIt]->getProcessor());
                        }
                    }
                }
            }
        }

        //if processor in actual layer was also in a lower layer, erase it from lower layer
        for (size_t procIt = 0; procIt < connProc.size(); procIt++) {
            for (size_t layerIt = 0; layerIt < layerList.size(); layerIt++) {
                std::vector<Processor*> lastLayer = layerList[layerIt];
                for (size_t it = 0; it < lastLayer.size(); it++) {
                    if (connProc[procIt] == lastLayer[it]) {
                        lastLayer.erase(lastLayer.begin() + it);
                        layerList.erase(layerList.begin() + layerIt);
                        layerList.insert(layerList.begin() + layerIt, lastLayer);
                    }
                }
            }
        }

        //add next layer to layerList
        layerList.push_back(connProc);

        //if there are no more processors to sort stop sorting and delete last (empty) entry
        if (connProc.empty()) {
            lastConn = true;
            layerList.pop_back();
        }

        //reset helping variables
        connProc.clear();
        ports.clear();
    }
    return layerList;

    //alternative sorting via GraphNodes
    //IN DEV!
    //
    /*
    std::vector<NetworkGraph::GraphNode*> nodes = netGraph.getGraphNodes();
    int numInQueue = 0;
    int newNumInQueue = 0;
    std::queue<Processor*> queue;
    std::vector<Processor*> procQueue;
    std::vector<std::vector<Processor*>> levelList;
    std::vector<Processor*> innerList;
    std::vector<NetworkGraph::GraphNode*> succNodes;
    std::queue<NetworkGraph::GraphNode*> save;

    for (size_t i = 0; i < nodes.size(); i++) {
        if (nodes[i]->getInDegree() == 0) {
            save.push(nodes[i]);
            queue.push(nodes[i]->getProcessor());
            procQueue.push_back(nodes[i]->getProcessor());
            numInQueue++;
        }
    }
    int count = 0;
    while (numInQueue) {
        for (int i = 0; i < numInQueue; i++) {
            //for (int j = 0; j < queue.size(); j++) {
                innerList.push_back(queue.front());
                queue.pop();
        }
        levelList.push_back(innerList);

        for (size_t i = 0; i < save.size(); i++) {
            succNodes = save.front()->getSuccessors();
            for (size_t j = 0; j < succNodes.size(); j++) {
                int in = succNodes[j]->getInDegree()-1;
                if(in == 0){
                    queue.push(succNodes[j]->getProcessor());
                }
            }
            save.pop();
        }
        for (size_t i = 0; i < succNodes.size(); i++) {
            save.push(succNodes[i]);
        }
        succNodes.clear();

        newNumInQueue = queue.size();
        numInQueue = newNumInQueue;
        newNumInQueue = 0;

    }
    return levelList;
    */
}

std::vector<Processor*> NWEGL_Base::getRenderingOrder(ProcessorNetwork* network) const {

    // construct graph: ignore loop ports for graph construction (second parameter), but unroll them afterwards (third parameter)
    NetworkGraph netGraph(network->getProcessors(), LoopPortTypeCheck(true), LoopPortTypeCheck());

    // generate subset of processors that are predecessors of the end processors
    std::vector<Processor*> endProcessors = getEndProcessors(network);
    std::set<Processor*> predecessors = netGraph.getPredecessors(endProcessors);

    // sort generated subset topologically
    std::vector<Processor*> renderingOrder = netGraph.sortTopologically(predecessors);

    return renderingOrder;
}

std::vector<Processor*> NWEGL_Base::getConnectedToOut(Processor* proc) {

    std::vector<Processor*> connProc;

    std::vector<Port*> out = proc->getOutports();
    for (size_t i = 0; i < out.size(); i++) {
        std::vector<const Port*> in = out[i]->getConnected();
        for (size_t j = 0; j < in.size(); j++) {
            Processor* connect = in[j]->getProcessor();
            connProc.push_back(connect);
        }
    }
    return connProc;
}

std::vector<Processor*> NWEGL_Base::getConnectedToIn(Processor* proc) {

    std::vector<Processor*> connProc;

    std::vector<Port*> in = proc->getInports();
    for (size_t i = 0; i < in.size(); i++) {
        std::vector<const Port*> out = in[i]->getConnected();
        for (size_t j = 0; j < out.size(); j++) {
            Processor* connect = out[j]->getProcessor();
            connProc.push_back(connect);
        }
    }
    return connProc;
}

size_t NWEGL_Base::getLayer(Processor* proc, std::vector<std::vector<Processor*> > graphLayer) {

    //if graphLayer is empty return -1
    if (graphLayer.empty())
        return -1;

    else {
        for(size_t layerIt = 0; layerIt < graphLayer.size(); layerIt++) {
            std::vector<Processor*> actLayer = graphLayer[layerIt];
                if(std::find(actLayer.begin(), actLayer.end(), proc)!=actLayer.end()) {
                    return layerIt;
                }
                //if processor is not found in graphLayer return -1
                else return -1;
        }
    }
    return -1;
}

} // namespace voreen
