/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_CONNECTEDCOMPONENTS3D_H
#define VRN_CONNECTEDCOMPONENTS3D_H

#include "voreen/core/processors/volumeprocessor.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"

namespace voreen {

/**
 * Detects connected components in a volume data set and intensity-codes the assigned labels
 * in a 16 bit output volume of the same dimensions.
 *
 * @see ConnectedComponents2D
 */
class ConnectedComponents3D : public CachingVolumeProcessor {
public:
    ConnectedComponents3D();
    ~ConnectedComponents3D();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "ConnectedComponents3D"; }
    virtual std::string getCategory() const     { return "Volume Processing";     }
    virtual CodeState getCodeState() const      { return CODE_STATE_STABLE;       }

protected:
    virtual void setDescriptions() {
        setDescription("Detects connected components in a volume data set and intensity-codes the assigned labels in a 16 bit output volume of the same dimensions.\
<p><strong>Properties</strong>:\
<ul>\
 <li>Connectivity: voxel neighborhood to consider.</li>\
 <li>Min Component Size: components consisting of less voxels are discarded.</li>\
 <li>Max Components: maximal number of connected components to put out.</li>\
 <li>Component Sorting: determines in which order the component labels are assigned.</li>\
 <li>Binarize Output: all connected components are assigned the same label (max intensity).</li>\
 <li>Stretch Labels: component labels are stretched over the whole intensity range of the output volume in order to maximize contrast.</li>\
</ul></p>\
<p>See: ConnectedComponents2D</p>");
    }

    virtual void process();

private:
    VolumePort inport_;     ///< Volume to analyze.
    VolumePort outport_;    ///< Output volume storing the components' labels.

    BoolProperty enableProcessing_;         ///< If set to false, the input volume is passed through.
    IntOptionProperty connectivity_;        ///< Voxel neighborhood to consider for the analysis.
    IntProperty minComponentSize_;          ///< Components consisting of less voxels are discarded.
    IntProperty maxComponents_;             ///< The maximal number of connected components that are put out.
    StringOptionProperty componentSorting_; ///< Determines in which order the component labels are assigned.
    BoolProperty binarizeOutput_;           ///< If enabled, all connected components are assigned the same label (max intensity).
    BoolProperty stretchLabels_;            ///< If enabled, component labels are stretched over the whole intensity range
                                            ///  of the output volume in order to maximize contrast.

    static const std::string loggerCat_;
};

} // namespace

#endif // VRN_CONNECTEDCOMPONENTS3D_H
