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

    vtkSmartPointer<vtkXMLMultiBlockDataReader> reader = vtkXMLMultiBlockDataReader::New();
    if(reader->CanReadFile(urlOrigin.getPath().c_str())) {
        reader->SetFileName(urlOrigin.getPath().c_str());
        reader->Update();

        vtkSmartPointer<vtkMultiBlockDataSet> blockData = vtkMultiBlockDataSet::SafeDownCast(reader->GetOutput());

        int blockIdx = 0;
        vtkSmartPointer<vtkCompositeDataIterator> iter = blockData->NewIterator();
        for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextItem()) {
            vtkSmartPointer<vtkImageData> block = vtkImageData::SafeDownCast(iter->GetCurrentDataObject());

            for (int i = 0; i < block->GetPointData()->GetNumberOfArrays(); i++) {
                VolumeURL subURL("vtm", url, "");
                const char *channel = block->GetPointData()->GetArrayName(i);
                subURL.addSearchParameter("scalar", channel);
                subURL.addSearchParameter("block", std::to_string(blockIdx));
                result.push_back(subURL);
            }

            blockIdx++;
        }
    }

    return result;
}

VolumeBase* VTMVolumeReader::read(const VolumeURL& origin) {

    std::string fileName = origin.getPath();
    LINFO("Reading " << origin.getURL());

    vtkSmartPointer<vtkXMLMultiBlockDataReader> reader = vtkXMLMultiBlockDataReader::New();

    if(!reader->CanReadFile(fileName.c_str())) {
        throw tgt::IOException();
    }

    reader->SetFileName(fileName.c_str());
    reader->Update();

    vtkSmartPointer<vtkMultiBlockDataSet> blockData = vtkMultiBlockDataSet::SafeDownCast(reader->GetOutput());
    unsigned int blockIdx = static_cast<unsigned int>(std::stoi(origin.getSearchParameter("block")));
    tgtAssert(blockIdx < blockData->GetNumberOfBlocks(), "Invalid block index");

    vtkSmartPointer<vtkMultiBlockDataSet> block = vtkMultiBlockDataSet::SafeDownCast(blockData->GetBlock(blockIdx));
    if(!block || block->GetNumberOfBlocks() == 0)
        throw tgt::IOException("Block could not be read.");

    if(block->GetNumberOfBlocks() > 1)
        LWARNING("Dataset contains more than one block. Ignoring.");

    return createVolumeFromVtkImageData(origin, vtkImageData::SafeDownCast(block->GetBlock(0)));
}

VolumeList* VTMVolumeReader::read(const std::string &url) {
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

VolumeReader* VTMVolumeReader::create(ProgressBar* progress) const {
    return new VTMVolumeReader(progress);
}

} // namespace voreen
