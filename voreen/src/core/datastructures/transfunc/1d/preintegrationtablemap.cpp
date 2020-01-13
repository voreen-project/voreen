/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "voreen/core/datastructures/transfunc/1d/preintegrationtablemap.h"

#include "voreen/core/datastructures/transfunc/1d/preintegrationtable.h"
#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"

namespace voreen {

PreIntegrationTableMap::PreIntegrationTableMap(TransFunc1D* tf, size_t maxSize) : tf_(tf) {
    //maxSize_ must be >= 1
    if (maxSize == 0)
        maxSize_ = 1;
    else
        maxSize_ = maxSize;
}

PreIntegrationTableMap::~PreIntegrationTableMap() {
    clear();
}

void PreIntegrationTableMap::setTransFunc(TransFunc1D* tf) {
    tf_ = tf;
}

void PreIntegrationTableMap::clear() {
    //delete all pre-integration tables, then clear map and queue
    std::map<PreIntegrationTableAttributes, PreIntegrationTable*>::iterator it = tableMap_.begin();
    for (; it!=tableMap_.end(); ++it)
        delete it->second;

    tableMap_.clear();
    queue_ = std::queue<PreIntegrationTableAttributes>();
}

const PreIntegrationTable* PreIntegrationTableMap::getPreIntegrationTable(float samplingStepSize, size_t dimension, bool useIntegral, bool computeOnGPU, tgt::Shader* program) {

    if (!tf_)
        return 0;

    size_t dim;

    if (dimension == 0)
        dim = static_cast<size_t>(std::min(tf_->getDimensions().x, 1024));
    else
        dim = dimension;

    //create key to pre-integration table map
    PreIntegrationTableAttributes attributes;
    attributes.samplingStepSize_ = samplingStepSize;
    attributes.dimension_ = dim;
    attributes.useIntegral_ = useIntegral;
    attributes.computeOnGPU_ = computeOnGPU;

    //try to find the right pre-integration table
    std::map<PreIntegrationTableAttributes, PreIntegrationTable*>::iterator it = tableMap_.find(attributes);
    if (it != tableMap_.end())
        return (it->second);
    else {
        //if size == max: delete the first
        if (tableMap_.size() == maxSize_) {
            it = tableMap_.find(queue_.front());
            delete it->second;
            tableMap_.erase(it);
            queue_.pop();
        }
        //create new pre-integration table and add it to the map
        PreIntegrationTable* table = new PreIntegrationTable(tf_, dim, samplingStepSize, useIntegral, computeOnGPU, program);
        tableMap_.insert(std::make_pair(attributes,table));
        queue_.push(attributes);
        it = tableMap_.find(attributes);
        return (it->second);
    }
}

} //namespace
