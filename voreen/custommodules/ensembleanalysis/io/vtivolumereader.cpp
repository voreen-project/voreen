/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#include "vtivolumereader.h"

#include <limits>

#include <vtkDataArray.h>
#include <vtkAbstractArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataReader.h>

#include "tgt/exception.h"
#include "tgt/logmanager.h"
#include "tgt/vector.h"

#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"

namespace voreen {

const std::string VTIVolumeReader::loggerCat_ = "voreen.io.VolumeReader.vti";

VTIVolumeReader::VTIVolumeReader(ProgressBar* progress)
    : VolumeReader(progress)
{
    extensions_.push_back("vti");
    protocols_.push_back("vti");
}

std::vector<VolumeURL> VTIVolumeReader::listVolumes(const std::string& url) const {
    std::vector<VolumeURL> result;

    VolumeURL urlOrigin(url);

    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkXMLImageDataReader::New();
    if(reader->CanReadFile(urlOrigin.getPath().c_str())) {
        reader->SetFileName(urlOrigin.getPath().c_str());
        reader->Update();

        for(int i = 0; i < reader->GetNumberOfPointArrays(); i++) {
            VolumeURL subURL("vti", url, "");
            const char* channel = reader->GetPointArrayName(i);
            subURL.addSearchParameter("scalar", channel);
            result.push_back(subURL);
        }
    }

    return result;
}

VolumeBase* VTIVolumeReader::read(const VolumeURL& origin) {

    std::string fileName = origin.getPath();
    LINFO("Reading " << origin.getURL());

    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkXMLImageDataReader::New();
    if(!reader->CanReadFile(fileName.c_str())) {
        throw tgt::IOException();
    }

    reader->SetFileName(fileName.c_str());
    reader->Update();

    vtkSmartPointer<vtkImageData> imageData = reader->GetOutput();
    tgt::svec3 dimensions = tgt::ivec3::fromPointer(imageData->GetDimensions());
    tgt::vec3 spacing = tgt::dvec3::fromPointer(imageData->GetSpacing());
    tgt::vec3 offset = tgt::dvec3::fromPointer(imageData->GetOrigin());
    std::string channel = origin.getSearchParameter("scalar");
    float min, max;

    VolumeRAM* dataset;
    try {
        // Try to retrieve data array.
        vtkDataArray* array = imageData->GetPointData()->GetArray(channel.c_str());
        if(!array)
            throw tgt::IOException("Field " + channel + " could not be read.");

        // Allocate volume.
        switch (array->GetNumberOfComponents()) {
        case 1:
            dataset = new VolumeRAM_Float(dimensions);
            break;
        case 2:
            dataset = new VolumeRAM_2xFloat(dimensions);
            break;
        case 3:
            dataset = new VolumeRAM_3xFloat(dimensions);
            break;
        case 4:
            dataset = new VolumeRAM_4xFloat(dimensions);
            break;
        default:
            throw tgt::IOException("Unsupported number of components");
        }

        // Export data.
        array->ExportToVoidPointer(dataset->getData());

        // Extract range.
        double range[2];
        array->GetRange(range);
        min = static_cast<float>(range[0]);
        max = static_cast<float>(range[1]);

    } catch (std::bad_alloc&) {
        throw; // throw it to the caller
    }

    Volume* volumeHandle = new Volume(dataset, spacing, offset);
    if (volumeHandle->getNumChannels() == 1) {
        volumeHandle->addDerivedData(new VolumeMinMax(min, max, min, max));
    }
    else {
        volumeHandle->addDerivedData(new VolumeMinMaxMagnitude(min, max));
    }
    volumeHandle->setOrigin(origin);
    volumeHandle->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping(volumeHandle->getBaseType()));
    volumeHandle->getMetaDataContainer().addMetaData("Scalar", new StringMetaData(channel));

    // Read meta data.
    for(int i = 0; i < imageData->GetFieldData()->GetNumberOfArrays(); i++) {
        vtkAbstractArray* array = imageData->GetFieldData()->GetAbstractArray(i);

        MetaDataBase* metaData = nullptr;
        switch (array->GetDataType()) {
        case VTK_INT:
            metaData = new IntMetaData(vtkIntArray::FastDownCast(array)->GetValue(0));
            break;
        case VTK_FLOAT:
            metaData = new FloatMetaData(vtkFloatArray::FastDownCast(array)->GetValue(0));
            break;
        case VTK_DOUBLE:
            metaData = new DoubleMetaData(vtkDoubleArray::FastDownCast(array)->GetValue(0));
            break;
        default:
            //LWARNING("Unsupported Meta Data found: " << array->GetName());
            break;
        }

        if(!metaData)
            continue;

        volumeHandle->getMetaDataContainer().addMetaData(array->GetName(), metaData);
    }

    return volumeHandle;
}

VolumeList* VTIVolumeReader::read(const std::string &url) {
    std::vector<VolumeURL> urls = listVolumes(url);
    VolumeList* volumeList = new VolumeList();

    try {
        for(const VolumeURL& url : urls) {
            volumeList->add(read(url));
        }
    } catch(tgt::FileException e) {
        while(!volumeList->empty()) {
            delete volumeList->first();
            volumeList->remove(volumeList->first());
        }
        delete volumeList;
        throw;
    }
    return volumeList;
}

VolumeReader* VTIVolumeReader::create(ProgressBar* progress) const {
    return new VTIVolumeReader(progress);
}

} // namespace voreen
