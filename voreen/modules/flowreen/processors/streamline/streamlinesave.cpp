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

#include "streamlinesave.h"

#include "voreen/core/voreenapplication.h"
#include "modules/core/io/vvdformat.h"

#include "tgt/filesystem.h"

#include <fstream>

namespace voreen {

const std::string StreamlineSave::loggerCat_("voreen.flowreen.StreamlineSave");

StreamlineSave::StreamlineSave()
    : Processor()
    // ports
    , streamlineInport_(Port::INPORT, "streamlineInport", "Streamline Input", false)
    , volumeInport_(Port::INPORT, "volumeInport", "Volume input", false)
    // properties
    , filenameProp_("filenameprop", "Save File as", "Select file...", VoreenApplication::app()->getUserDataPath(),
            "Voreen streamline data (*.vsd);;Character-separated values (*.csv)", FileDialogProperty::SAVE_FILE, Processor::INVALID_PATH)
    , saveButton_("saveButton", "Save")
    , saveVolumeInVsd_("saveVolumeInVsd","Save Volume in VSD?",true,Processor::INVALID_RESULT,Property::LOD_ADVANCED)
    // members
    , saveStreamlines_(false)
{
    addPort(streamlineInport_);
    addPort(volumeInport_);

    saveButton_.onChange(MemberFunctionCallback<StreamlineSave>(this, &StreamlineSave::saveStreamlines));

    addProperty(filenameProp_);
    addProperty(saveButton_);
    addProperty(saveVolumeInVsd_);
}

StreamlineSave::~StreamlineSave() {
}

bool StreamlineSave::isReady() const {
    if (!isInitialized())
        return false;
    // only streamlineport must be ready
    if(!streamlineInport_.isReady())
            return false;

    return true;
}

void StreamlineSave::invalidate(int inv) {
    Processor::invalidate(inv);

    if(inv == Processor::INVALID_PATH && isInitialized()) {
        saveStreamlines_ = true;
    }
}

void StreamlineSave::process() {
    if (saveStreamlines_){
        saveStreamlines();
        saveStreamlines_ = false;
    }
}

//---------------------------------------------------------------------------------
//      Callbacks
//---------------------------------------------------------------------------------
void StreamlineSave::saveStreamlines() {
    if (!isInitialized())
        return;
    if (!streamlineInport_.getData()) {
        LWARNING("no input streamlines");
        return;
    }
    if (filenameProp_.get().empty()) {
        LWARNING("no filename specified");
        return;
    }
    std::string extension = FileSys.fileExtension(filenameProp_.get(),true);
    if(extension.compare("csv") && extension.compare("vsd")) {
        LWARNING("Selected file extension is neither *.csv nor *.vsd");
        return;
    }
    // write file
    try {
        std::ofstream outFile;
        outFile.open(filenameProp_.get().c_str());
        LINFO("Writing streamlines to file " << filenameProp_.get());

        if(!extension.compare("csv")) {
            //write meta
            outFile << streamlineInport_.getData()->metaToCSVString().c_str() << "\n";
            //write streamlines
            for(size_t i = 0; i < streamlineInport_.getData()->getStreamlines().size(); i++) {
                outFile << streamlineInport_.getData()->getStreamlines().at(i).toCSVString(streamlineInport_.getData()->getVoxelToWorldMatrix(),
                                                                                           streamlineInport_.getData()->getVelocityTransformMatrix()).c_str() << "\n";
            }
        } else { // case vsd
            XmlSerializer s;
            Serializer serializer(s);
            serializer.serialize("VSD-File", streamlineInport_.getData());
            bool saveVolume = false;
            //save volume if present
            if(saveVolumeInVsd_.get()) {
                if(!volumeInport_.hasData()) {
                    LWARNING("No volume in \"Volume Input\", although saveVolumeInVsd is selected. Saving without volume.");
                } else {
                    saveVolume = true;
                }
            }
            if(saveVolume) {
                VvdObject vvd(volumeInport_.getData(),"",false);
                std::vector<VvdObject> vec;
                vec.push_back(vvd);
                serializer.serialize("Magnitude-Volumes", vec, "Volume");
                s.write(outFile); // vvd must present during write. null pointer otherwise
            } else {
                s.write(outFile);
            }
        }

        outFile.close();
    }
    catch(tgt::FileException e) {
        LERROR(e.what());
        filenameProp_.set("");
    }
    LINFO("Saving " << filenameProp_.get() << " was successful.");
}

}   // namespace
