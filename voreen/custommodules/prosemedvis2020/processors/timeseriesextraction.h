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

#ifndef VRN_TIMESERIESEXTRACTION_H
#define VRN_TIMESERIESEXTRACTION_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/string/stringlistproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"

#include "../modules/ensembleanalysis/ports/ensembledatasetport.h"
#include "../modules/plotting/ports/plotport.h"
#include "../modules/plotting/datastructures/plotdata.h"
#include "../ports/ensemblesamplepointmappingport.h"

#include <vector>

namespace voreen {
    
/**
 *  Processor, which extracts a timeseries of PET-Data from a given point.
 */
class VRN_CORE_API TimeseriesExtraction : public Processor {
public:

    TimeseriesExtraction();
    virtual Processor* create() const;
    virtual std::string getClassName() const;
    virtual std::string getCategory() const;

protected:

    virtual void setDescriptions();
    virtual void process();

private:

    // Ports
    EnsembleDatasetPort ensembleInport_;                    ///< inport used to receive the ensemble
    EnsembleSamplePointMappingPort samplePointInport_;      ///< inport used to recieve the sample-points
    PlotPort outport_;                                      ///< outport used to pass the plot-data

    // Propertys
    OptionProperty<std::string> fieldProp_;// Property used to select a field from the ensemble
    BoolProperty useAgvSampling_;          // Used to enable and disable avg-sampling
    FloatProperty sphereDiameter_;         // Diameter of the sphere used for avg-sampling
    OptionProperty<float> timeUnitProp_; // Property to select the output time-unit
    StringProperty selTimeUnitProp_;    // Displays the selected Time-unit (uset to link with plot axis label)
    OptionProperty<std::string> unitProp_; // Property to select the output unit
    StringProperty selUnitProp_;    // Displays the selected unit (uset to link with plot axis label)

    /** Called when the input of ensembleInport_ changes */
    void ensembleChanged();

    /** Returns an averaged sample-value unsing a sphere of the diameter specified in sphereDiameter_  
     *  The spacing-vector is used to match the size of the sphere to the size of rendered Sample-Points.
    */
    float getAvgSample(const tgt::vec3& samplepoint, const VolumeRAMRepresentationLock& volumeData, const tgt::vec3& spacing);

    float convertTime(float timeInSeconds);

    const float getTotalCount(const VolumeRAMRepresentationLock& volumeData, const RealWorldMapping& rwm, const tgt::vec3& spacing) const;
};

} // namespace

#endif // VRN_TIMESERIESEXTRACTION_H
