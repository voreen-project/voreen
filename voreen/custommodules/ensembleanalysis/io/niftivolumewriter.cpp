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

#include "niftivolumewriter.h"

#include <vtkAbstractArray.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkNIFTIImageWriter.h>

#include "tgt/exception.h"
#include "tgt/vector.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

const std::string NiftiVolumeWriter::loggerCat_ = "voreen.io.VolumeWriter.nii";

NiftiVolumeWriter::NiftiVolumeWriter(ProgressBar* progress)
        : VolumeWriter(progress)
{
    extensions_.push_back("nii");
}

VolumeWriter* NiftiVolumeWriter::create(ProgressBar* progress) const {
    return new NiftiVolumeWriter(progress);
}

void NiftiVolumeWriter::write(const std::string& fileName, const VolumeBase* volumeHandle) {
    tgtAssert(volumeHandle, "No volume");

    LINFO("Writing " << fileName);

    vtkSmartPointer<vtkNIFTIImageWriter> writer = vtkSmartPointer<vtkNIFTIImageWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(createVtkImageDataFromVolume(volumeHandle));
    writer->Write();
    if(writer->GetErrorCode() != 0) {
        throw tgt::IOException("File could not be written");
    }
}

} // namespace voreen
