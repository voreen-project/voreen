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

#ifndef VRN_SEGMENTATIONTHRESHOLDTOTRANSFUNCLINKER_H
#define VRN_SEGMENTATIONTHRESHOLDTOTRANSFUNCLINKER_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/volumeport.h"

#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

namespace voreen {

    /**
     * This class is used to link the octreegrowing processor with a sliceviewer.
     * The threshold defined by the octreegrowing processor will be visualized in
     * the sliceviewer.
     */
class VRN_CORE_API ThresholdToTransFuncLinker : public Processor {
public:
    ThresholdToTransFuncLinker();
    //virtual ~ThresholdToTransFuncLinker();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "ThresholdToTransFuncLinker";      }
    virtual std::string getCategory() const     { return "Linker";  }
    virtual CodeState getCodeState() const      { return CODE_STATE_EXPERIMENTAL;  }

    virtual bool isReady() const;

protected:
    /**
     * Enum used to specify which threshold value should be vsualized.
     */
    enum StttflMode {
        STTTFL_QUANTIFICATION_A,
        STTTFL_QUANTIFICATION_B,
        STTTFL_QUANTIFICATION_C
    };

    virtual void initialize();

    //virtual void beforeProcess();
    virtual void process();
    //virtual void afterProcess();

    virtual void setDescriptions() {
        setDescription(" This class is used to link the octreegrowing processor with a sliceviewer. \
                         The threshold defined by the octreegrowing processor will be visualized in \
                         the sliceviewer.");
    }

private:
    /**
     * The function sets the tf of property \"channel\" to the color of the associated color
     * property. The tf must be of class TransFunc1DKeys.
     * @param channel The tf property to modify. Values between 1-4
     */
    void onChangeCreateOpaqueTransFunc(size_t channel);

    /**
     * Function used to adjust all tf properties to the current channel and threshold.
     */
    void onChangeThresholdChannel();

    /**
     * If the parameter is the current value of the thresholdChannel_ property onChangeThresholdChannel()
     * will be called.
     */
    void onChangeOtherProperty(StttflMode associatedMode);

    void adjustToInputVolume();

    VolumePort volumePort_; ///< port used to get the real world mapping for the tf properties

    OptionProperty<StttflMode> thresholdChannel_; ///< determines the currently visualized channel

    TransFunc1DKeysProperty channel1TransFuncProp_;   ///< linked tf to sliceviewer
    ColorProperty color1Prop_;                  ///< color of the tf
    TransFunc1DKeysProperty channel2TransFuncProp_;   ///< linked tf to sliceviewer
    ColorProperty color2Prop_;                  ///< color of the tf
    TransFunc1DKeysProperty channel3TransFuncProp_;   ///< linked tf to sliceviewer
    ColorProperty color3Prop_;                  ///< color of the tf
    TransFunc1DKeysProperty channel4TransFuncProp_;   ///< linked tf to sliceviewer
    ColorProperty color4Prop_;                  ///< color of the tf

    IntProperty channelAProp_;                      ///< linked prop to octreequantification
    IntIntervalProperty thresholdAProp_;            ///< linked prop to octreequantification
    IntProperty maxValueAProp_;                     ///< linked prop to octreequantification

    IntProperty channelBProp_;                      ///< linked prop to octreequantification
    IntIntervalProperty thresholdBProp_;            ///< linked prop to octreequantification
    IntProperty maxValueBProp_;                     ///< linked prop to octreequantification

    IntProperty channelCProp_;                      ///< linked prop to octreequantification
    IntIntervalProperty thresholdCProp_;            ///< linked prop to octreequantification
    IntProperty maxValueCProp_;                     ///< linked prop to octreequantification

    static const std::string loggerCat_;        ///< miau
};

} // namespace

#endif // VRN_THRESHOLDTOTRANSFUNCLINKER_H
