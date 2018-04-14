/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_NWEGL_BASE_H
#define VRN_NWEGL_BASE_H

#include "voreen/qt/voreenqtapi.h"

#include <vector>
#include <queue>

//qt
#include <QMap>

namespace voreen {
    class ProcessorNetwork;
    class Processor;
    class ProcessorGraphicsItem;
    class NetworkGraph;
    class GraphNode;

class VRN_QT_API NWEGL_Base {
public:
    NWEGL_Base() {}
    virtual ~NWEGL_Base() {}

    //function which sorts the graph
    virtual void sort(ProcessorNetwork* network, std::vector<Processor*>* processors, QMap<Processor*,ProcessorGraphicsItem*>* itemMap) = 0;

protected:
    //functions from the evaluator
    std::vector<Processor*> getSourceProcessors(ProcessorNetwork* network) const;
    std::vector<Processor*> getEndProcessors(ProcessorNetwork* network) const;
    std::vector<Processor*> getRenderingOrder(ProcessorNetwork* network) const;

    //assigns every processor of the network to a layer (vertical sorted) and return
    std::vector<std::vector<Processor*> > getGraphLayer(std::vector<Processor*> proc) const;

    //returns a vector of processors that are connected to an outport of a processor
    std::vector<Processor*> getConnectedToOut(Processor* proc);

    //returns a vector of processors that are connected to an inport of a processor
    std::vector<Processor*> getConnectedToIn(Processor* proc);

    //get Layer of processor in given Layerlist
    size_t getLayer(Processor* proc, std::vector<std::vector<Processor*> > graphLayer);
};

} //namespace voreen

#endif // VRN_NWEGL_BASE_H

