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

#ifndef VRN_STREAMLINESAVE_H
#define VRN_STREAMLINESAVE_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/boolproperty.h"

#include "../../ports/streamlinelistport.h"
#include "voreen/core/ports/volumeport.h"

namespace voreen {

    /**
     * Processor used to save streamlines calculated from the StreamlineCreator to Disk.
     * Supported file formats are CSV and Voreen Streamline Data (VSD).
     * If a Volume is connected as inport, the volume will be stored in the VSD file.
     *
     * @see StreamlineCreator
     * @see Streamline
     */
class VRN_CORE_API StreamlineSave : public Processor {
public:
    /** Constructior */
    StreamlineSave();
    /** Destructor */
    virtual ~StreamlineSave();

    virtual Processor* create() const         { return new StreamlineSave(); }
    virtual std::string getClassName() const  { return "StreamlineSave";     }
    virtual std::string getCategory() const   { return "Output";             }
    virtual CodeState getCodeState() const    { return CODE_STATE_STABLE;    }
    virtual bool isEndProcessor() const       { return true;                 }

protected:
    virtual void setDescriptions() {
        setDescription("Processor used to save streamlines calculated from the StreamlineCreator to Disk. " \
                       "Supported file formats are CSV and Voreen Streamline Data (VSD). The CSV formatation is:<br>" \
                       "Dimension, Bounds.LLF, Bounds.URB<br>Streamline1<br>Streamline2<br>...<br>" \
                       "Bundle1<br>Bundle2<br>...<br>" \
                       "Where each streamline is stored as:<br>position1, velocity1, ... positionX, velocityX<br><br>" \
                       "and each bundle is stored as:<br>radius, position1, velocity1, ... positionX, velocityX<br><br>" \
                       "The VSD format can be loaded into Voreen by the StreamlineSource. If the volume port is connected, the " \
                       "volume will be stored in the vsd file for later usage.");
        streamlineInport_.setDescription("Port containing the streamline data to be saved.");
        volumeInport_.setDescription("If a volume is connected and the save volume property is enabled,the volume will be stored in the vsd file. "\
                                     "For a csv file, this port has no functionality. Normally this file corresponds the the magnitude image used as reference.");
    }

    //--------------------------
    //  Override
    //--------------------------
    /** Saves the volume file. */
    virtual void process() override;
    /** Processor should be ready, even if no volume is attached */
    virtual bool isReady() const override;
    /** Triggers a auto-save if the file property has changed values. */
    virtual void invalidate(int inv = 1) override;

    /** Supported file formats */
    enum StreamlineFormat {
        SF_CSV,
        SF_VSD
    };

    //--------------
    //  Callbacks
    //--------------
    /** Main function used to save the streamlines. Triggered by saveButton_. */
    void saveStreamlines();

    //--------------
    //  Member
    //--------------
private:
    StreamlineListPort streamlineInport_;       ///< inport containing the streamlines to save
    VolumePort volumeInport_;                   ///< magnitude volume, which can be stored in the vsd file

    FileDialogProperty filenameProp_;               ///< determines the name of the saved file
    ButtonProperty saveButton_;                     ///< triggers a save
    BoolProperty saveVolumeInVsd_;                  ///< if true, a connected volume is stored inside the vsd file

    bool saveStreamlines_;          ///< used to determine, if process should save or not

    static const std::string loggerCat_;
};

}   //namespace

#endif
