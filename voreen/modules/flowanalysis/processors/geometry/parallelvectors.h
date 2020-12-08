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

#ifndef VRN_PARALLELVECTORS_H
#define VRN_PARALLELVECTORS_H

#include <string>
#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"

#include "modules/flowanalysis/ports/parallelvectorsolutionpointsport.h"

namespace voreen {

/**
 * This processor implements the parallel vectors operator by Peikert and Roth and optional sujudi-haimes filtering.
 */
class ParallelVectors : public Processor {
public:
    ParallelVectors();
    virtual Processor* create() const { return new ParallelVectors(); }
    virtual std::string getClassName() const { return "ParallelVectors"; }
    virtual std::string getCategory() const { return "Volume Processing"; }
    virtual bool isReady() const;

    static void Process(const VolumeRAM& v, const VolumeRAM& w, const VolumeRAM_Mat3Float* jacobian, const VolumeRAM* mask, ParallelVectorSolutions& outSolution, const RealWorldMapping& rwmV = RealWorldMapping(), const RealWorldMapping& rwmW = RealWorldMapping());

    static constexpr auto TetrahedraPerCube = 6;
    static constexpr auto TrianglesPerTetrahedron = 4;

protected:
    virtual void process();

    void setDescriptions() {
        setDescription("This processor implements the parallel vectors operator by Peikert and Roth and optional sujudi-haimes filtering. "
                       "It can be used to extract vortex corelines by using velocity, acceleration volumes, as well as the jacobi matrix as input.");
        _inV.setDescription("First input volume (V) for the parallel vectors operator");
        _inW.setDescription("Second input volume (W) for the parallel vectors operator");
        _inJacobian.setDescription("(Optional) Jacobi matrix to be used for sujudi-haimes filtering");
    }

private:
    void onChangedJacobianData();

    VolumePort _inV, _inW, _inJacobian, _inMask;
    ParallelVectorSolutionPointsPort _out;
    BoolProperty _sujudiHaimes;

    using Triangle = std::array<tgt::svec3, 3>;
    using Tet = std::array<Triangle, TrianglesPerTetrahedron>;
};

} // namespace voreen

#endif