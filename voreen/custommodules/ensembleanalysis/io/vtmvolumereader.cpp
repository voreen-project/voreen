/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "vtmvolumereader.h"

#include "vtivolumereader.h"

#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkCompositeDataIterator.h>

#include <vtkSmartPointer.h>
#include <vtkXMLMultiBlockDataReader.h>
#include <vtkInformation.h>

#include "tgt/exception.h"
#include "tgt/logmanager.h"

namespace voreen {

const std::string VTMVolumeReader::loggerCat_ = "voreen.io.VolumeReader.vtm";

VTMVolumeReader::VTMVolumeReader(ProgressBar* progress)
    : VolumeReader(progress)
{
    extensions_.push_back("vtm");
    protocols_.push_back("vtm");
}

std::vector<VolumeURL> VTMVolumeReader::listVolumes(const std::string& url) const {
    std::vector<VolumeURL> result;

    VolumeURL urlOrigin(url);
    std::string filterName = urlOrigin.getSearchParameter("name");
    std::string filterBlock = urlOrigin.getSearchParameter("block");

    vtkSmartPointer<vtkXMLMultiBlockDataReader> reader = vtkSmartPointer<vtkXMLMultiBlockDataReader>::New();
    if(reader->CanReadFile(urlOrigin.getPath().c_str())) {
        reader->SetFileName(urlOrigin.getPath().c_str());
        reader->Update();

        vtkMultiBlockDataSet* blockData = vtkMultiBlockDataSet::SafeDownCast(reader->GetOutput());

        int blockIdx = 0;
        vtkCompositeDataIterator* iter = blockData->NewIterator();
        for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextItem()) {
            vtkImageData* block = vtkImageData::SafeDownCast(iter->GetCurrentDataObject());

            std::string blockIdxString = std::to_string(blockIdx);
            if(!filterBlock.empty() && filterBlock != blockIdxString)
                continue;

            for (int i = 0; i < block->GetPointData()->GetNumberOfArrays(); i++) {

                const char* name = block->GetPointData()->GetArrayName(i);
                if (!filterName.empty() && filterName != name)
                    continue;

                VolumeURL subURL("vtm", urlOrigin.getPath(), "");
                subURL.addSearchParameter("name", name);
                subURL.addSearchParameter("block", blockIdxString);
                subURL.getMetaDataContainer().addMetaData("name", new StringMetaData(name));
                result.push_back(subURL);
            }

            blockIdx++;
        }
        iter->Delete();
    }

    return result;
}

VolumeBase* VTMVolumeReader::read(const VolumeURL& origin) {

    std::string fileName = origin.getPath();
    LINFO("Reading " << origin.getURL());

    vtkSmartPointer<vtkXMLMultiBlockDataReader> reader = vtkSmartPointer<vtkXMLMultiBlockDataReader>::New();
    if(!reader->CanReadFile(fileName.c_str())) {
        throw tgt::IOException("File can't be read!");
    }

    reader->SetFileName(fileName.c_str());
    reader->Update();

    vtkMultiBlockDataSet* blockData = vtkMultiBlockDataSet::SafeDownCast(reader->GetOutput());
    unsigned int blockIdx = static_cast<unsigned int>(std::stoi(origin.getSearchParameter("block")));
    tgtAssert(blockIdx < blockData->GetNumberOfBlocks(), "Invalid block index");

    vtkMultiBlockDataSet* block = vtkMultiBlockDataSet::SafeDownCast(blockData->GetBlock(blockIdx));
    if(!block || block->GetNumberOfBlocks() == 0)
        throw tgt::IOException("Block could not be read.");

    if(block->GetNumberOfBlocks() > 1)
        LWARNING("Dataset contains more than one block. Ignoring.");

    return createVolumeFromVtkImageData(origin, vtkImageData::SafeDownCast(block->GetBlock(0)));
}

VolumeList* VTMVolumeReader::read(const std::string& url) {
    std::vector<VolumeURL> urls = listVolumes(url);
    std::vector<std::unique_ptr<VolumeBase>> volumes;

    try {
        for(const auto& url : urls) {
            volumes.push_back(std::unique_ptr<VolumeBase>(read(url)));
        }
    } catch(tgt::IOException& e) {
        throw;
    }

    // Transfer ownership to output list.
    VolumeList* volumeList = new VolumeList();
    for(auto& volume : volumes) {
        volumeList->add(volume.release());
    }
    return volumeList;
}

VolumeReader* VTMVolumeReader::create(ProgressBar* progress) const {
    return new VTMVolumeReader(progress);
}

} // namespace voreen
