/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_AORTASEGMENTATION_H
#define VRN_AORTASEGMENTATION_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

#include "voreen/core/ports/volumeport.h"

namespace voreen {

class HeartSegmentationBorder {
public:
    HeartSegmentationBorder(const tgt::vec3& aortaPos, const tgt::vec3& heartPos, float radius);

    bool isCrossedBy(const tgt::ivec3& before, const tgt::ivec3& after) const;

private:
    bool isInSphereOfInfluence(const tgt::vec3&) const;
    bool areOnOppositeSidesOfPlane(const tgt::vec3& p1, const tgt::vec3& p2) const;

    tgt::vec3 center_;
    tgt::vec3 planeNormal_;
    float effectRadiusSquared_;
    float cylinderHeight_;
};

/**
 * Interactive volume segmentation by region growing.
 *
 * Input:
 *  - volume to be segmented
 *  - geometry with two points: First has to be on/in the vessel, second has to be in/on the heart
 * Output: segmented volume
 */
class VRN_CORE_API AortaSegmentation : public Processor {
public:
    AortaSegmentation();

    virtual ~AortaSegmentation();

    virtual std::string getCategory() const { return "Volume Processing"; }
    virtual std::string getClassName() const { return "AortaSegmentation"; }
    virtual Processor* create() const { return new AortaSegmentation(); }

    virtual void process();
    virtual void adjustPropertiesToInput();

protected:
    virtual void setDescriptions() {
        setDescription("Volume segmentation by region growing.");
    }

    void mark(const tgt::ivec3& seedpos, int segment, const HeartSegmentationBorder& border, VolumeRAM_UInt8& segVol, const VolumeBase& volumeInput, const VolumeBase& vesselnessInput);

    FloatProperty strictness_;         ///< influences cost function for filling
    BoolProperty thresholdFilling_;    ///< use threshold as additional restriction for filling?
    FloatIntervalProperty volumeThresholds_;     ///< intensity thresholds which are applied when thresholdFilling_ is true
    FloatIntervalProperty vesselnessThresholds_;     ///< intensity thresholds which are applied when thresholdFilling_ is true
    StringOptionProperty fillCostFunction_;    ///< cost function to use
    BoolProperty adaptive_;            ///< adjust criteria during growing
    FloatProperty separatingDiskRadius_;
    IntProperty maxSteps_;
    BoolProperty updateContinuously_;            ///< adjust criteria during growing
    ButtonProperty startGrowing_;

    bool updateForced_;

    GeometryPort seedPoints_;
    GeometryPort heartPoints_;
    VolumePort volumeInport_;
    VolumePort vesselnessInport_;
    VolumePort outport_;

    static const std::string loggerCat_; ///< category used in logging
};

} // namespace

#endif // VRN_AORTASEGMENTATION_H
