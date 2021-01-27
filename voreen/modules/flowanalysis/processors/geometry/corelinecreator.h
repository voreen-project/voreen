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

#ifndef VRN_CORELINECREATOR_H
#define VRN_CORELINECREATOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/intproperty.h"

#include "modules/flowanalysis/ports/parallelvectorsolutionpointsport.h"

namespace voreen {

class CorelineCreator : public Processor {
public:
    CorelineCreator();
    virtual Processor* create() const { return new CorelineCreator(); }
    virtual std::string getClassName() const { return "CorelineCreator"; }
    virtual std::string getCategory() const { return "Geometry"; }
    virtual CodeState getCodeState() const { return CODE_STATE_TESTING; }

    static void Process( const ParallelVectorSolutions& solutions, int lengthThreshold, std::vector<std::vector<tgt::vec3>>& corelines );

protected:
    virtual void process();

private:
    ParallelVectorSolutionPointsPort _in;
    GeometryPort _out;
    IntProperty _lengthThreshold;
};

} // namespace voreen

#endif