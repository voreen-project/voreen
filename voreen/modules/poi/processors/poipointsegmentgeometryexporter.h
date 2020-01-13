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

#ifndef VRN_POIPOINTSEGMENTEXPORTER
#define VRN_POIPOINTSEGMENTEXPORTER

#include "voreen/core/processors/processor.h"
#include "poistorage.h"
#include "voreen/core/ports/genericcoprocessorport.h"
#include "voreen/core/ports/geometryport.h"
#include <voreen/core/datastructures/geometry/pointsegmentlistgeometry.h>
namespace voreen {

/**
 * Converts Points of interest to a PointSegmentGeometryVec3, with each group as a segment
 */
class VRN_CORE_API POIPointSegmentExporter : public Processor {

public:
    POIPointSegmentExporter();

    virtual Processor* create() const;
    virtual std::string getClassName() const;
    virtual std::string getCategory() const;
protected:
    virtual void setDescriptions();
    virtual void process();
private:
    GenericCoProcessorPort<POIStorage> cpPort_;
    GeometryPort geometryport_;
    std::unique_ptr<PointSegmentListGeometryVec3> geometry_;
};

} // namespace

#endif
