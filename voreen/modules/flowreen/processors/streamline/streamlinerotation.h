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

#ifndef VRN_STREAMLINEROTATION_H
#define VRN_STREAMLINEROTATION_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/boolproperty.h"

#include "voreen/core/datastructures/volume/slice/slicetexture.h"

#include "../../ports/streamlinelistport.h"
#include "voreen/core/ports/volumeport.h"

namespace voreen {

    /**
     * Processor used to rotate and align volumes/streamlines according to an PlaneManipulation
     */
class VRN_CORE_API StreamlineRotation : public Processor {
public:
    StreamlineRotation();
    virtual ~StreamlineRotation();
    virtual Processor* create() const override        { return new StreamlineRotation(); }
    virtual std::string getClassName() const override { return "StreamlineRotation";  }
    virtual std::string getCategory() const override  { return "Process";          }
    virtual CodeState getCodeState() const override   { return CODE_STATE_STABLE; }

protected:
    virtual void setDescriptions() override {
        setDescription("This Processor can be used to rotate and align Streamlines. For a proper use, a StreamlineCreater for the streamlines, a \
                        PlaneManipulation for the configuration and a VolumeTransform for the volume transformation are required.");
    }

    virtual void process() override;
    virtual void serialize(Serializer& s) const override;
    virtual void deserialize(Deserializer& d) override;

    //--------------
    //  Callbacks
    //--------------
    /** StreamlinePort gets new data*/
    void onNewDataStreamlinePort();

    /** Called, if the rotation axis is updated. */
    void rotationAxisOnChange();
    /** Called, if the rotation angle is updated. */
    void rotationAngleOnChange();
    /** Called, if the linked handle is updated. */
    void linkedHandleOnChange();


    /** Updates the linked matrix. */
    void updateLinkedMatrix();

    /** Updates the bool flags. */
    void alignmentSettingsOnChange();
    /** Adjusts the alignmentMatrix_ */
    void applyAlignmentOnChange();
    /** Clears the alignmentMatrix_ */
    void clearAlignmentOnChange();

    //--------------
    //  Member
    //--------------
private:
    //ports
    StreamlineListPort streamlineInport_;   ///< inport used in process
    StreamlineListPort streamlineOutport_;  ///< outport for the rotated lines

    //rotation
    FloatVec3Property rotationAxisProp_;    ///< axis used for the rotation
    FloatProperty     rotationAngleProp_;   ///< angle used for the rotation
    BoolProperty linkAxisAndHandleProp_;    ///< if true, the handle and the rotation axis are linked
    //alignment
    BoolProperty alignmentIsAppliedProp_;                   ///< alignment in general is applied
    BoolProperty currentAlignmentSettingsAreAppliedProp_;   ///< shows, if the current settings are applied
    FloatVec3Property linkedHandleNormalProp_;              ///< this normal will be aligned
    OptionProperty<SliceAlignment> alignToPlaneProp_;       ///< aligned to this plane
    ButtonProperty applyAlignmentProp_;                     ///< apply alignment
    ButtonProperty clearAlignmentProp_;                     ///< clears the current alignment

    tgt::mat4 transMatOfLastStreamlineList_;            ///< stores the transformation matrix of the last input streamline ist
    tgt::mat4 currentRotationMatrix_;                   ///< current used roation. A combination of rotationAxisProp_ and rotationAngleProp_
    tgt::mat4 alignmentMatrix_;                         ///< current matrix of an rotation from the handle onto the alignment plane normal
    tgt::mat4 lastConfigurationVelocityRotation_;       ///< last used configuration of rotation/alignment for the Streamlines. Only the roation part
    tgt::mat4 lastConfigurationStreamlineTransform_;    ///< last used configuration of rotation/alignment for the streamline list

    //linked result
    FloatMat4Property linkedMatrixProp_;    ///< rotation matrix

    static const std::string loggerCat_;    ///< meow
};

}   //namespace

#endif
