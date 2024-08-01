/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_PREINTEGRATIONTABLEMAP_H
#define VRN_PREINTEGRATIONTABLEMAP_H

#include "voreen/core/voreencoreapi.h"

#include "tgt/vector.h"
#include "tgt/shadermanager.h"

#include <queue>
#include <map>

namespace voreen {

class PreIntegrationTable;
class TransFunc1D;

/// Used as key value for the PreIntegrationTableMap class.
struct PreIntegrationTableAttributes {
    float samplingStepSize_;
    size_t dimension_;
    bool useIntegral_;
    bool computeOnGPU_;

    // operator< needed for inserting as keys into std::map
    bool operator< (const PreIntegrationTableAttributes& a) const {
        if (samplingStepSize_ < a.samplingStepSize_)
            return true;
        else if (samplingStepSize_ == a.samplingStepSize_) {
            if (dimension_ < a.dimension_)
                return true;
            else if (dimension_ == a.dimension_) {
                if (useIntegral_ && !a.useIntegral_)
                    return true;
                else if (useIntegral_ == a.useIntegral_) {
                    if (computeOnGPU_ && !a.computeOnGPU_)
                        return true;
                }
            }
        }

        return false;
    }
};

/**
 * A map that holds several pre-integration tables for the same transfer function.
 */
class VRN_CORE_API PreIntegrationTableMap {

public:

    PreIntegrationTableMap(TransFunc1D* tf, size_t maxSize = 1);

    ~PreIntegrationTableMap();

    void setTransFunc(TransFunc1D* tf);

    /**
     * Returns the pre-integration table with the specified attributes.
     * If necessary, the table is (re-)computed.
     * If the table map already contains the maximum number of tables, the one inserted first is deleted from the table (FIFO).
     *
     * @param dimension width of the pre-integration table, 0 chooses the width according to the bit depth of the volume (up to 1024)
     * @param samplingStepSize the segment length used for rendering
     * @param useIntegral @see PreIntegrationTable
     */
    const PreIntegrationTable* getPreIntegrationTable(float samplingStepSize = 1.f, size_t dimension = 0, bool useIntegral = true, bool computeOnGPU = false, tgt::Shader* program = 0);

    /**
     * Delete alle pre-integration tables in the map (e.g. if tf has changed).
     */
    void clear();

private:

    PreIntegrationTableMap();

    TransFunc1D* tf_; ///< transfer function to compute the pre-integration tables

    std::map<PreIntegrationTableAttributes, PreIntegrationTable*> tableMap_; ///< map containing the pre-integration tables
    std::queue<PreIntegrationTableAttributes> queue_; ///< stores attributes of the tables in insertion order

    size_t maxSize_; ///< maximum number of elements allowed in the map

};

} //namespace

#endif // VRN_PREINTEGRATIONTABLEMAP_H
