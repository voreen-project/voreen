/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "netcdfvolumereader.h"

#include "vtivolumereader.h"

#include <vtkCellData.h>
#include <vtkImageData.h>
#include <vtkNetCDFCFReader.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>

#include "tgt/exception.h"
#include "tgt/logmanager.h"

namespace voreen {

const std::string NetCDFVolumeReader::loggerCat_ = "voreen.io.VolumeReader.netcdf";

NetCDFVolumeReader::NetCDFVolumeReader(ProgressBar* progress)
    : VolumeReader(progress)
{
    extensions_.push_back("nc");
    protocols_.push_back("nc");
}

std::vector<VolumeURL> NetCDFVolumeReader::listVolumes(const std::string& url) const {
    std::vector<VolumeURL> result;

    VolumeURL urlOrigin(url);

    vtkSmartPointer<vtkNetCDFCFReader> reader = vtkSmartPointer<vtkNetCDFCFReader>::New();
    if(reader->CanReadFile(urlOrigin.getPath().c_str())) {
        reader->SetFileName(urlOrigin.getPath().c_str());
        reader->UpdateMetaData();

        int numArrays = reader->GetNumberOfVariableArrays();
        for(int i=0; i<numArrays; i++) {
            const char* name = reader->GetVariableArrayName(i);
            VolumeURL subURL("nc", urlOrigin.getPath(), "");
            subURL.addSearchParameter("name", name);
            //subURL.addSearchParameter("timeStep", std::to_string(t)); // Seems not to be possible using NetCDF CF.
            subURL.getMetaDataContainer().addMetaData("name", new StringMetaData(name));
            result.push_back(subURL);
        }
    }
    else {
        throw tgt::IOException("File seems not to follow the CF convention");
    }

    return result;
}

VolumeBase* NetCDFVolumeReader::read(const VolumeURL& origin) {

    std::string fileName = origin.getPath();
    LINFO("Reading " << origin.getURL());

    vtkSmartPointer<vtkNetCDFCFReader> reader = vtkSmartPointer<vtkNetCDFCFReader>::New();
    if(!reader->CanReadFile(fileName.c_str())) {
        throw tgt::IOException("File can't be read!");
    }

    reader->SetFileName(fileName.c_str());
    reader->UpdateMetaData();

    vtkStringArray* variableDimensions = reader->GetVariableDimensions();
    for(vtkIdType i=0; i<variableDimensions->GetNumberOfValues(); i++) {
        reader->SetVariableArrayStatus(reader->GetVariableArrayName(i), 0);
    }

    std::string name = origin.getSearchParameter("name");
    reader->SetVariableArrayStatus(name.c_str(), 1);
    //double timeStep = std::stod(origin.getSearchParameter("timeStep")); // Seems not to be possible using NetCDF CF.
    //reader->UpdateTimeStep(timeStep);
    reader->SetOutputTypeToImage();
    reader->Update();

    vtkImageData* imageData = vtkImageData::SafeDownCast(reader->GetOutput());
    //imageData->GetFieldData()->AddArray(); // Add meta data here.
    Volume* volume = createVolumeFromVtkImageData(origin, imageData);

    // Query unit, if available.
    RealWorldMapping rwm = volume->getRealWorldMapping();
    rwm.setUnit(reader->QueryArrayUnits(name.c_str()));
    volume->setRealWorldMapping(rwm);

    return volume;
}

VolumeList* NetCDFVolumeReader::read(const std::string& url) {
    std::vector<VolumeURL> urls = listVolumes(url);
    std::vector<std::unique_ptr<VolumeBase>> volumes;

    try {
        for(const auto& url : urls) {
            volumes.push_back(std::unique_ptr<VolumeBase>(read(url)));
        }
    } catch(tgt::IOException&) {
        throw;
    }

    // Transfer ownership to output list.
    VolumeList* volumeList = new VolumeList();
    for(auto& volume : volumes) {
        volumeList->add(volume.release());
    }
    return volumeList;
}

VolumeReader* NetCDFVolumeReader::create(ProgressBar* progress) const {
    return new NetCDFVolumeReader(progress);
}

bool NetCDFVolumeReader::canSupportFileWatching() const {
    return true;
}

} // namespace voreen
