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

#include "vtivolumewriter.h"

#include <vtkAbstractArray.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>

#include "tgt/exception.h"
#include "tgt/vector.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

vtkSmartPointer<vtkImageData> createVtkImageDataFromVolume(const VolumeBase* volume) {
    VolumeRAMRepresentationLock representation(volume);

    // Setup image data object.
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    tgt::ivec3 dims = representation->getDimensions();
    imageData->SetDimensions(dims.x, dims.y, dims.z);
    tgt::dvec3 spacing = volume->getSpacing();
    imageData->SetSpacing(spacing.x, spacing.y, spacing.z);
    tgt::dvec3 origin = volume->getOffset();
    imageData->SetOrigin(origin.x, origin.y, origin.z);
    int numComponents = static_cast<int>(representation->getNumChannels());

    int dataType = VTK_VOID;
    if(representation->getBaseType() == "int8") {
        dataType = VTK_CHAR;
    }
    else if(representation->getBaseType() == "uint8") {
        dataType = VTK_UNSIGNED_CHAR;
    }
    else if(representation->getBaseType() == "int16") {
        dataType = VTK_SHORT;
    }
    else if(representation->getBaseType() == "uint16") {
        dataType = VTK_UNSIGNED_SHORT;
    }
    else if(representation->getBaseType() == "int32") {
        dataType = VTK_INT;
    }
    else if(representation->getBaseType() == "uint32") {
        dataType = VTK_UNSIGNED_INT;
    }
    else if(representation->getBaseType() == "int64") {
        dataType = VTK_LONG;
    }
    else if(representation->getBaseType() == "uint64") {
        dataType = VTK_UNSIGNED_LONG;
    }
    else if(representation->getBaseType() == "float") {
        dataType = VTK_FLOAT;
    }
    else if(representation->getBaseType() == "double") {
        dataType = VTK_DOUBLE;
    }
    else {
        throw tgt::UnsupportedFormatException("VTK format not supported.");
    }
    imageData->AllocateScalars(dataType, numComponents);

    RealWorldMapping rwm = volume->getRealWorldMapping();
    for(int z=0; z<dims.z; z++) {
        for(int y=0; y<dims.y; y++) {
            for(int x=0; x<dims.x; x++) {
                for(int component = 0; component < numComponents; component++) {
                    float value = rwm.normalizedToRealWorld(representation->getVoxelNormalized(x, y, z, component));
                    imageData->SetScalarComponentFromFloat(x, y, z, component, value);
                }
            }
        }
    }

    return imageData;
}


const std::string VTIVolumeWriter::loggerCat_ = "voreen.io.VolumeWriter.vti";

VTIVolumeWriter::VTIVolumeWriter(ProgressBar* progress)
        : VolumeWriter(progress)
{
    extensions_.push_back("vti");
}

VolumeWriter* VTIVolumeWriter::create(ProgressBar* progress) const {
    return new VTIVolumeWriter(progress);
}

void VTIVolumeWriter::write(const std::string& fileName, const VolumeBase* volumeHandle) {
    tgtAssert(volumeHandle, "No volume");

    LINFO("Writing " << fileName);

    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(createVtkImageDataFromVolume(volumeHandle));
    if(!writer->Write()) {
        throw tgt::IOException("File could not be written");
    }
}

} // namespace voreen
