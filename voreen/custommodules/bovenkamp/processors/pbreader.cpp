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

#include "pbreader.h"
#include "../io/volumediskpb.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/filesystem.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumelist.h"

#include <iostream>
#include <fstream>
#include <sstream>

#ifdef UNIX
#include <clocale>
#endif

namespace voreen {

const std::string PBReader::loggerCat_("voreen.bovenkamp.PBReader");

std::string FILE_PARAMETER = "FOV_parameters.txt";
std::string FILE_MAGNITUDE = "magnitude.txt";
std::string FILE_VELOCITYX = "velocity_x.txt";
std::string FILE_VELOCITYY = "velocity_y.txt";
std::string FILE_VELOCITYZ = "velocity_z.txt";

PBReader::PBReader()
    : Processor()
    , magnitudeOutport_(Port::OUTPORT,"magintudeOutport","Magnitude")
    , velocityOutport_(Port::OUTPORT,"velocityOutport","Velocitye")
    , folderProp_("folderProp","Data Folder","Select Folder", VoreenApplication::app()->getUserDataPath(),"",FileDialogProperty::DIRECTORY,Processor::VALID)
    , loadButtonProp_("loadButtonProp", "Load Files",Processor::INVALID_RESULT)
    , invertXInputProp_("invertXInputProp","Invert Input X",false,Processor::VALID)
    , invertYInputProp_("invertYInputProp","Invert Input Y",true,Processor::VALID)
    , invertZInputProp_("invertZInputProp","Invert Input Z",true,Processor::VALID)
    , invertXVelocityProp_("invertXVelocityProp","Invert Velocity X",false,Processor::VALID)
    , invertYVelocityProp_("invertYVelocityProp","Invert Velocity Y",false,Processor::VALID)
    , invertZVelocityProp_("invertZVelocityProp","Invert Velocity Z",false,Processor::VALID)
    , loadFiles_(false)
    , isMagnitudeDataPresent_(false)
    , isVelocityDataPresent_(false)
    {
    // add ports
    addPort(magnitudeOutport_);
    addPort(velocityOutport_);

    // add properties
    addProperty(folderProp_);
    addProperty(loadButtonProp_);
    loadButtonProp_.onChange(MemberFunctionCallback<PBReader>(this,&PBReader::load));

    addProperty(invertXInputProp_);
    addProperty(invertYInputProp_);
    addProperty(invertZInputProp_);
    addProperty(invertXVelocityProp_);
    addProperty(invertYVelocityProp_);
    addProperty(invertZVelocityProp_);

    invertXInputProp_.setGroupID("input");
    invertYInputProp_.setGroupID("input");
    invertZInputProp_.setGroupID("input");
    setPropertyGroupGuiName("input","Invert Input");
    invertXVelocityProp_.setGroupID("velocity");
    invertYVelocityProp_.setGroupID("velocity");
    invertZVelocityProp_.setGroupID("velocity");
    setPropertyGroupGuiName("velocity","Invert Velocity Direction");

    invertXInputProp_.onChange(MemberFunctionOneParameterCallback<PBReader,pbDirection>(this,&PBReader::invertInputOnChange,PB_X));
    invertYInputProp_.onChange(MemberFunctionOneParameterCallback<PBReader,pbDirection>(this,&PBReader::invertInputOnChange,PB_Y));
    invertZInputProp_.onChange(MemberFunctionOneParameterCallback<PBReader,pbDirection>(this,&PBReader::invertInputOnChange,PB_Z));
    invertXVelocityProp_.onChange(MemberFunctionCallback<PBReader>(this,&PBReader::invertVelocityOnChange));
    invertYVelocityProp_.onChange(MemberFunctionCallback<PBReader>(this,&PBReader::invertVelocityOnChange));
    invertZVelocityProp_.onChange(MemberFunctionCallback<PBReader>(this,&PBReader::invertVelocityOnChange));

}

PBReader::~PBReader() {
}

void PBReader::initialize() {
    Processor::initialize();
    //trigger the loading
    load();
}

bool PBReader::isReady() const {
    if(magnitudeOutport_.isConnected() || velocityOutport_.isConnected())
        return true;

    return false;
}

void PBReader::process() {
    if(loadFiles_) {

        //reset everything
        loadFiles_ = false;
        magnitudeVolumes_.clear();
        magnitudeOutport_.clear();
        velocityVolumes_.clear();
        velocityOutport_.clear();

        //get parameters
        tgt::svec3 dimensions;
        tgt::vec3 spacing;
        int timesteps;
        try {
            readParameters(dimensions, spacing, timesteps);
        }
        catch(VoreenException& e) {
            LERROR(e.what());
            return;
        }

        //load magnitude
        tgt::bvec3 invertPosition(invertXInputProp_.get(), invertYInputProp_.get(), invertZInputProp_.get());
        VolumeList* magnitudeList = new VolumeList();
        if(isMagnitudeDataPresent_) {
            for (int t = 0; t < timesteps; t++) {
                VolumeDiskPB* magnitude = new VolumeDiskPB(folderProp_.get() + "/" + FILE_MAGNITUDE, invertPosition, dimensions, t);
                Volume* volume = new Volume(magnitude, spacing, tgt::vec3::zero);
                magnitudeList->add(volume);
                magnitudeVolumes_.push_back(std::unique_ptr<Volume>(volume));
            }
        }
        magnitudeOutport_.setData(magnitudeList);

        //load velocity
        VolumeList* velocityList = new VolumeList();
        if(isVelocityDataPresent_) {
            tgt::bvec3 invertVelocity(invertXVelocityProp_.get(), invertYVelocityProp_.get(), invertZVelocityProp_.get());
            for (int t = 0; t < timesteps; t++) {
                VolumeDiskPB* velocity = new VolumeDiskPB(
                    folderProp_.get() + "/" + FILE_VELOCITYX,
                    folderProp_.get() + "/" + FILE_VELOCITYY,
                    folderProp_.get() + "/" + FILE_VELOCITYZ,
                    invertPosition, invertVelocity,
                    dimensions, t);
                Volume* volume = new Volume(velocity, spacing, tgt::vec3::zero);
                velocityList->add(volume);
                velocityVolumes_.push_back(std::unique_ptr<Volume>(volume));
            }
        }
        velocityOutport_.setData(velocityList);
    }
}

//---------------------------------------------------
//      Callbacks
//---------------------------------------------------
void PBReader::invertInputOnChange(pbDirection dir) {
    switch(dir) {
    case PB_X:
        invertXVelocityProp_.set(!invertXVelocityProp_.get());
        break;
    case PB_Y:
        invertYVelocityProp_.set(!invertYVelocityProp_.get());
        break;
    case PB_Z:
        invertZVelocityProp_.set(!invertZVelocityProp_.get());
        break;
    }
    //clearOutput is called in invertVelocityOnChange()
}

void PBReader::invertVelocityOnChange() {
    clearOutput();
}

void PBReader::clearOutput() {
    //clear evrything
    loadFiles_ = false;
    magnitudeVolumes_.clear();
    magnitudeOutport_.clear();
    velocityVolumes_.clear();
    velocityOutport_.clear();
    isMagnitudeDataPresent_ = false;
    isVelocityDataPresent_ = false;
}

//---------------------------------------------------
//      The Loading
//---------------------------------------------------
void PBReader::load() {
    isMagnitudeDataPresent_ = false;
    isVelocityDataPresent_ = false;

    //no folder => return
    if(folderProp_.get().empty()) {
        LWARNING("No data folder is selected!");
        return;
    }

    if (!tgt::FileSystem::dirExists(folderProp_.get())) {
        LWARNING("Directory does not exists: " << folderProp_.get());
        return;
    }

    LINFO("Loading from directory " << folderProp_.get());

    //no param => return
    if(!tgt::FileSystem::fileExists(folderProp_.get() + "/" + FILE_PARAMETER)) {
        LWARNING("Selected folder " << folderProp_.get() << " does not contain " << FILE_PARAMETER << "!");
        return;
    }

    //check if magnitude is present
    if(!tgt::FileSystem::fileExists(folderProp_.get() + "/" + FILE_MAGNITUDE)) {
        LWARNING("Selected folder " << folderProp_.get() << " does not contain " << FILE_MAGNITUDE << "! No magnitude volume will be loaded.");
    } else {
        isMagnitudeDataPresent_ = true;
    }

    //check if velocity is present
    isVelocityDataPresent_ = true;
    if(!tgt::FileSystem::fileExists(folderProp_.get() + "/" + FILE_VELOCITYX)) {
        LWARNING("Selected folder " << folderProp_.get() << " does not contain " << FILE_VELOCITYX << "! No velocity volume will be loaded.");
        isVelocityDataPresent_ = false;
    } else if(!tgt::FileSystem::fileExists(folderProp_.get() + "/" + FILE_VELOCITYY)) {
        LWARNING("Selected folder " << folderProp_.get() << " does not contain " << FILE_VELOCITYY << "! No velocity volume will be loaded.");
        isVelocityDataPresent_ = false;
    } else if(!tgt::FileSystem::fileExists(folderProp_.get() + "/" + FILE_VELOCITYZ)) {
        LWARNING("Selected folder " << folderProp_.get() << " does not contain " << FILE_VELOCITYZ << "! No velocity volume will be loaded.");
        isVelocityDataPresent_ = false;
    }

    //return if magnitude and velocity is missing
    if(!isMagnitudeDataPresent_ && !isVelocityDataPresent_) return;

    //file will be loaded during next process
    loadFiles_ = true;
}

void PBReader::readParameters(tgt::svec3& dimensions, tgt::vec3& spacing, int& timesteps) {
    //open file
    std::ifstream ifs((folderProp_.get() + "/" + FILE_PARAMETER).c_str());
    if (ifs.fail()) {
        throw VoreenException("Parameter file: '" + folderProp_.get() + "/" + FILE_PARAMETER + "' could not be opened");
    }
    //TODO: checks
    std::string tmp;
    getline(ifs, tmp); //number of pixels
    getline(ifs, tmp); //dimensions
    std::stringstream line(tmp);
    for(int d = 0; d < 3; d++) {
        getline(line,tmp,'\t');
        dimensions[d] = (size_t)atoi(tmp.c_str());
    }

    getline(ifs, tmp); //field of view
    getline(ifs, tmp); //values
    getline(ifs, tmp); //step size
    getline(ifs, tmp); //spacing
    line.str(tmp);
    for(int s = 0; s < 3; s++) {
        getline(line,tmp,'\t');
        spacing[s] = (float)atof(tmp.c_str());
    }
    getline(ifs, tmp); //number of time frames
    getline(ifs, tmp); //timesteps
    line.str(tmp);
    getline(line,tmp,'\t');
    timesteps = atoi(tmp.c_str());

    ifs.close();
};

} // namespace voreen





