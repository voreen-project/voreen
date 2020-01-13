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

#ifndef VRN_RADARGLYPHVOLUMEGENERATOR_H
#define VRN_RADARGLYPHVOLUMEGENERATOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/optionproperty.h"
namespace voreen {

class CameraInteractionHandler;

class RadarGlyphVolumeGenerator : public Processor
{
public:
    RadarGlyphVolumeGenerator();
    virtual ~RadarGlyphVolumeGenerator();
    virtual Processor* create() const { return new RadarGlyphVolumeGenerator(); }

    virtual std::string getCategory() const { return "Source"; }
    virtual std::string getClassName() const { return "RadarGlyphVolumeGenerator"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }
    virtual bool usesExpensiveComputation() const { return true; }
    virtual void process();
private:
    virtual void setDescriptions() {
        setDescription("");
    }

    enum GlyphOption {
        GO_2_3_4_20 = 0,
        GO_50_50_50_50 = 1,
        GO_50_50_50_500 = 2
    };
    friend class OptionProperty<int>;

    VolumeListPort outport_;

    bool generateNew_;
    OptionProperty<int> optionProp_;
    ButtonProperty generateProp_;
    void generateOnChange();
};

}   // namespace

#endif  // VRN_RADARGLYPHVOLUMEGENERATOR_H
