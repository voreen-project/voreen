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

#ifndef VRN_PBREADER_H
#define VRN_PBREADER_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/volumeport.h"

#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/filedialogproperty.h"

#include "voreen/core/utils/exception.h"

namespace voreen {

    /**
     * Processor for laoding velocity data provided by P. Bovenkamp.
     * TODO: Find clearTarget error in render port, if loading is canceled.
     */
class PBReader : public Processor {
public:
    PBReader();
    virtual ~PBReader();
    virtual Processor* create() const        { return new PBReader(); }
    virtual std::string getClassName() const { return "PBReader";  }
    virtual std::string getCategory() const  { return "Data Inport";     }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL; }
    virtual bool usesExpensiveComputation() const {return true;}

    /** Enables the laoding of the data for a new laoded workspace */
    void initialize();

protected:
    virtual void setDescriptions() {
        setDescription("Processor used to load files in the P.Bovenkamp format!");
    }

    virtual void process() override;
    virtual bool isReady() const override;

    /**
     * Checks, if the selected folder contains all needed files.
     * Sets loadFiles_ = true and stops running loadings.
     */
    void load();
    /**
     * Clears the Outports and stops running loadings.
     */
    void clearOutput();

    /** Enum to define invert direction. */
    enum pbDirection {
        PB_X = 0,
        PB_Y = 1,
        PB_Z = 2
    };

    void invertInputOnChange(pbDirection dir);
    void invertVelocityOnChange();
private:
    /**
     * Help function to read the parameters of the parameter file.
     * Called in "process".
     */
    void readParameters(tgt::svec3& dimensions, tgt::vec3& spacing, int& timesteps);
    /**
     * Help function to read the data of a magnitude or velocity file.
     * Called in "process".
     */
    void readFile(std::string filename, tgt::svec3 dimensions, int timesteps, float** outputArray, size_t currentChannel = 1, size_t numChannels = 3, 
                              float progressMin = 0.f, float progressMax = 1.f);
    
    /// Ports
    VolumeListPort magnitudeOutport_;   //< port containing the magnitude volume of all time steps
    VolumeListPort velocityOutport_;    //< port containing the velocity (vec3) volume of all time steps

    /// Properties
    FileDialogProperty folderProp_;     //< file dialog to select the folder containing the data files in pb-format
    ButtonProperty loadButtonProp_;     //< loads the data, if a valid folder is selected

    BoolProperty invertXInputProp_;      //< inverts the x input
    BoolProperty invertYInputProp_;      //< inverts the y input
    BoolProperty invertZInputProp_;      //< inverts the z input

    BoolProperty invertXVelocityProp_;      //< inverts the velocity direction for dimension x
    BoolProperty invertYVelocityProp_;      //< inverts the velocity direction for dimension y
    BoolProperty invertZVelocityProp_;      //< inverts the velocity direction for dimension z

    // Member
    bool loadFiles_;    //< if true, the process function will load the data set
    bool stopLoading_;  //< if true, the current loading is aborted

    bool isMagnitudeDataPresent_; //< load magnitude volume
    bool isVelocityDataPresent_;  //< load velocity volume

    static const std::string loggerCat_; //< static member for output LERROR etc...
};

} //namespace voreen

#endif
