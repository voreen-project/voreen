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

#ifndef VRN_TILEEXTRACTOR_H
#define VRN_TILEEXTRACTOR_H

#include "voreen/core/processors/volumeprocessor.h"

#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
//#include "voreen/core/properties/volumeinfoproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/stringproperty.h"

namespace voreen {

class VRN_CORE_API TileExtractor : public VolumeProcessor {

public:
    TileExtractor();

    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "TileExtractor";      }
    virtual std::string getCategory() const   { return "Classification";          }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }
    virtual bool isEndProcessor() const       { return true;              }

    virtual bool usesExpensiveComputation() const { return true; }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Exports 2D tiles with 3D information as channels as training data for Caffe in HDF5 format.");
    }

    virtual void process();

    virtual void computeTiles();

    virtual void deactivateProperties();

    virtual void reactivateProperties();

    tgt::vec2 getScalingMinMaxIntensity(const VolumeBase* volume) const;

    float scaleIntensity(float value, tgt::vec2 minMaxValue) const;

    /// computes maximum number of tiles in a single HDF5 file based on the 2GB limit set by Caffe 
    void adjustTilesPerFile();

private:

    // data input
    VolumeListPort dataInport_;

    // input of labels
    VolumeListPort labelInport_;
    
    // for scaling intensities
    VolumePort intensityReferenceVolume_;

    IntProperty tileSize_;
    IntProperty borderSize_;
    IntProperty channelWidth_;

    IntProperty numTilesPerFile_;

    BoolProperty useMirroring_;

    ButtonProperty computeButton_;

    StringProperty timeEstimation_;

    bool startComputation_;

    static const std::string loggerCat_;

    FileDialogProperty filename_;

    BoolProperty autoCompute_;

    ProgressProperty progressProperty_;
};

} // namespace
#endif 
