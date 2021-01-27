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

#ifndef VRN_PARALLELVECTORSOLUTIONPOINTSPORT_H
#define VRN_PARALLELVECTORSOLUTIONPOINTSPORT_H

#include "voreen/core/ports/genericport.h"
#include <vector>

namespace voreen {

struct ParallelVectorSolutions {
    std::vector<tgt::vec3> solutions; // Points at which vectors are parallel
    std::vector<int32_t> triangleSolutionIndices; // Indices of solutions inside each triangle or -1 if there is no solution
    tgt::svec3 dimensions; // Dimensions of the original volumes
    tgt::mat4 voxelToWorldMatrix; // Transformation matrix of the original volume.

    ParallelVectorSolutions() = default;
    ParallelVectorSolutions(std::vector<tgt::vec3>&& solutions, std::vector<int32_t>&& solutionIndices, tgt::svec3 dimensions, tgt::mat4 voxelToWorldMatrix)
        : solutions(std::move(solutions)), triangleSolutionIndices(std::move(solutionIndices)), dimensions(dimensions), voxelToWorldMatrix(voxelToWorldMatrix) {}
};

#ifdef DLL_TEMPLATE_INST
    template class VRN_CORE_API GenericPort<ParallelVectorSolutions>;
#endif

class VRN_CORE_API ParallelVectorSolutionPointsPort : public GenericPort<ParallelVectorSolutions> {
public:
    ParallelVectorSolutionPointsPort(PortDirection direction, const std::string& id, const std::string& guiName = "", bool allowMultipleConnections = false, Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT);

    virtual Port* create(PortDirection direction, const std::string& id, const std::string& guiName = "") const;
    virtual std::string getClassName() const;
    virtual tgt::col3 getColorHint() const;
    virtual std::string getContentDescription() const;
    virtual std::string getContentDescriptionHTML() const;
};

} // namespace voreen

#endif // VRN_PARALLELVECTORSOLUTIONPOINTSPORT_H
