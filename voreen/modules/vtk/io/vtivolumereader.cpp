/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "vtivolumereader.h"

#include "../vtkmodule.h"

#include <vtkAbstractArray.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataReader.h>

#include "tgt/exception.h"
#include "tgt/logmanager.h"
#include "tgt/vector.h"

#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"

#ifdef VRN_MODULE_HDF5
#include "modules/hdf5/io/hdf5volumereader.h"
#include "modules/hdf5/io/hdf5volumewriter.h"
#include "modules/hdf5/io/hdf5filevolume.h"
#endif

namespace voreen {

#ifdef VRN_MODULE_HDF5
class VolumeRAMSwap : public VolumeDisk {
public:

    static bool trySetVolumeSwapDisk(Volume* volume) {
        if(!volume) {
            return false;
        }

        std::unique_ptr<VolumeRAMSwap> swap;
        try {
            swap.reset(new VolumeRAMSwap(volume));
        } catch(tgt::IOException&) {
            return false;
        }

        volume->addRepresentation(swap.release());
        volume->removeRepresentation<VolumeRAM>();

        return true;
    }

    VolumeRAMSwap(Volume* volume)
        : VolumeDisk(volume->getFormat(), volume->getDimensions())
        , hash_(volume->getHash())
    {
        // Shorthands.
        const tgt::svec3 dim = volume->getDimensions();
        const std::string baseType = volume->getBaseType();
        const size_t numChannels = volume->getNumChannels();

        // File settings.
        const tgt::svec3 chunkSize(dim.xy(), 1);
        const int compressionLevel = 5;
        const bool shuffleEnabled = true;
        const std::string location = "swap";

        // Retrieve tmp name.
        std::string filename = VoreenApplication::app()->getUniqueTmpFilePath(".h5");

        // Write.
        fileVolume_ = std::unique_ptr<HDF5FileVolume>(HDF5FileVolume::createVolume(filename, location, baseType, dim, numChannels, true, compressionLevel, chunkSize, shuffleEnabled));
        fileVolume_->writeVolume(volume->getRepresentation<VolumeRAM>());
        fileVolume_->writeSpacing(volume->getSpacing());
        fileVolume_->writeOffset(volume->getOffset());
        fileVolume_->writePhysicalToWorldTransformation(volume->getPhysicalToWorldMatrix());
        if(volume->hasMetaData(VolumeBase::META_DATA_NAME_REAL_WORLD_MAPPING)) {
            fileVolume_->writeRealWorldMapping(volume->getRealWorldMapping());
        }
        fileVolume_->writeDerivedData(volume);
        fileVolume_->writeMetaData(volume);
    }

    std::string getHash() const {
        return hash_;
    }

    VolumeRAM* loadVolume() const {
        return fileVolume_->loadVolume(0, getNumChannels());
    }

    VolumeRAM* loadSlices(const size_t firstZSlice, const size_t lastZSlice) const {
        return fileVolume_->loadSlices(firstZSlice, lastZSlice, 0, getNumChannels());
    }

    VolumeRAM* loadBrick(const tgt::svec3& offset, const tgt::svec3& dimensions) const {
        return fileVolume_->loadBrick(offset, dimensions, 0, getNumChannels());
    }

private:

    std::unique_ptr<HDF5FileVolume> fileVolume_;
    std::string hash_;
};
#endif

Volume* createVolumeFromVtkImageData(const VolumeURL& origin, vtkSmartPointer<vtkImageData> imageData) {
    tgtAssert(imageData, "ImageData null");

    tgt::svec3 dimensions = tgt::ivec3::fromPointer(imageData->GetDimensions());
    tgt::vec3 spacing = tgt::dvec3::fromPointer(imageData->GetSpacing());
    tgt::vec3 offset = tgt::dvec3::fromPointer(imageData->GetOrigin());
    std::string name = origin.getSearchParameter("name");

    // Try to retrieve data array - favor cell data.
    vtkDataArray* array = imageData->GetCellData()->GetArray(name.c_str());
    if(!array) {
        array = imageData->GetPointData()->GetArray(name.c_str());
    }
    else {
        // FIXME: Somehow it seems we need to substract 1 in each dimension for cell data?
        dimensions -= tgt::svec3::one;
    }
    if(!array) {
        throw tgt::IOException("Field " + name + " could not be read.");
    }

    // Convert vtk data type to voreen data type.
    std::string dataType;
    switch (array->GetDataType()) {
    case VTK_CHAR:
    case VTK_SIGNED_CHAR:
        dataType = "int8";
        break;
    case VTK_UNSIGNED_CHAR:
        dataType = "uint8";
        break;
    case VTK_SHORT:
        dataType = "int16";
        break;
    case VTK_UNSIGNED_SHORT:
        dataType = "uint16";
        break;
    case VTK_INT:
        dataType = "int32";
        break;
    case VTK_UNSIGNED_INT:
        dataType = "uint32";
        break;
    case VTK_LONG:
        dataType = "int64";
        break;
    case VTK_UNSIGNED_LONG:
        dataType = "uint64";
        break;
    case VTK_FLOAT:
        dataType = "float";
        break;
    case VTK_DOUBLE:
        dataType = "double";
        break;
    default:
        throw tgt::UnsupportedFormatException("VTK format not supported.");
    }

    VolumeFactory volumeFactory;
    std::string format = volumeFactory.getFormat(dataType, array->GetNumberOfComponents());

    // Create volume.
    VolumeRAM* dataset = nullptr;
    try {
        dataset = volumeFactory.create(format, dimensions);
    } catch (std::bad_alloc&) {
        throw; // throw it to the caller
    }

    // Export data.
    array->ExportToVoidPointer(dataset->getData());

    // Extract range.
    std::vector<float> min;
    std::vector<float> max;
    for(int i=0; i<array->GetNumberOfComponents(); i++) {
        double range[2];
        array->GetRange(range, i);
        min.push_back(range[0]);
        max.push_back(range[1]);
    }

    Volume* volume = new Volume(dataset, spacing, offset);
    volume->addDerivedData(new VolumeMinMax(min, max, min, max));
    volume->setOrigin(origin);
    volume->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping(volume->getBaseType()));
    volume->setModality(Modality(name));

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

        volume->getMetaDataContainer().addMetaData(array->GetName(), metaData);
    }

#ifdef VRN_MODULE_HDF5
    // Add swap disk representation.
    if(VTKModule::getInstance() && VTKModule::getInstance()->getForceDiskRepresentation()) {
        VolumeRAMSwap::trySetVolumeSwapDisk(volume);
    }
#endif

    return volume;
}

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
    std::string filterName = urlOrigin.getSearchParameter("name");

    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
    if(reader->CanReadFile(urlOrigin.getPath().c_str())) {
        reader->SetFileName(urlOrigin.getPath().c_str());
        reader->UpdateInformation();

        for(int i = 0; i < reader->GetNumberOfCellArrays(); i++) {
            const char* name = reader->GetCellArrayName(i);

            if(!filterName.empty() && filterName != name)
                continue;

            VolumeURL subURL("vti", urlOrigin.getPath(), "");
            subURL.addSearchParameter("name", name);
            subURL.getMetaDataContainer().addMetaData("name", new StringMetaData(name));
            result.push_back(subURL);
        }

        for(int i = 0; i < reader->GetNumberOfPointArrays(); i++) {

            const char* name = reader->GetPointArrayName(i);

            if(!filterName.empty() && filterName != name)
                continue;

            VolumeURL subURL("vti", urlOrigin.getPath(), "");
            subURL.addSearchParameter("name", name);
            subURL.getMetaDataContainer().addMetaData("name", new StringMetaData(name));
            result.push_back(subURL);
        }
    }

    return result;
}

VolumeBase* VTIVolumeReader::read(const VolumeURL& origin) {

    std::string fileName = origin.getPath();
    LINFO("Reading " << origin.getURL());

    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
    if(!reader->CanReadFile(fileName.c_str())) {
        throw tgt::IOException("File can't be read!");
    }

    reader->SetFileName(fileName.c_str());
    reader->Update();

    return createVolumeFromVtkImageData(origin, reader->GetOutput());
}

VolumeList* VTIVolumeReader::read(const std::string& url) {
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

VolumeReader* VTIVolumeReader::create(ProgressBar* progress) const {
    return new VTIVolumeReader(progress);
}

bool VTIVolumeReader::canSupportFileWatching() const {
    return true;
}

} // namespace voreen
