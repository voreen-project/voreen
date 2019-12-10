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

#include "itkvolumereader.h"

#include <set>
#include <limits>

#include "tgt/exception.h"
#include "tgt/logmanager.h"
#include "tgt/vector.h"

#include "itkImageFileReader.h"
#include "itkImageIOBase.h"

#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "modules/itk/utils/itkwrapper.h"

using tgt::ivec3;
using tgt::vec3;

namespace voreen {

const std::string ITKVolumeReader::loggerCat_ = "voreen.io.VolumeReader.ITK";

template<class T>
    Volume* readScalarVolume(const std::string fileName) {
        typedef itk::Image<T,3> ImageNDType;
        typedef itk::ImageFileReader<ImageNDType> ReaderType;

        typename ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(fileName.c_str());
        int status = 1;
        try
        {
            reader->Update();
            status = 0;
        }
        catch(itk::ExceptionObject &excp)
        {
            throw tgt::FileException("ITK Exception: " + std::string(excp.GetDescription()), fileName);
        }
        if (status)
        {
            throw tgt::FileException("Other error reading file using ITK!", fileName);
        }

        return ITKToVoreenCopy<T>(reader->GetOutput());
    }

VolumeCollection* ITKVolumeReader::read(const std::string &url) {
    VolumeURL origin(url);
    std::string fileName = origin.getPath();

    LINFO("Reading file " << fileName);

    //Get OutputInformation of an arbitrary reader to find out pixel type etc:
    typedef itk::Image<char,3> TestImageType; // pixel type doesn't matter for current purpose
    typedef itk::ImageFileReader<TestImageType> TestFileReaderType; // reader for testing a file
    TestFileReaderType::Pointer onefileReader = TestFileReaderType::New();
    onefileReader->SetFileName(fileName.c_str());
    try
    {
        onefileReader->GenerateOutputInformation();
    }
    catch(itk::ExceptionObject& excp)
    {
        throw tgt::CorruptedFileException("Failed to read OutputInformation! " + std::string(excp.GetDescription()), fileName);
    }

    // grab the ImageIO instance for the reader
    itk::ImageIOBase *imageIO = onefileReader->GetImageIO();

    unsigned int NumberOfDimensions =  imageIO->GetNumberOfDimensions();
    LINFO("Number of Dimensions: " << NumberOfDimensions);
    if(NumberOfDimensions != 3) {
        throw tgt::UnsupportedFormatException("Unsupported number of dimensions!");
    }

    // PixelType is SCALAR, RGB, RGBA, VECTOR, COVARIANTVECTOR, POINT, INDEX
    itk::ImageIOBase::IOPixelType pixelType = imageIO->GetPixelType();
    LINFO("PixelType: " << imageIO->GetPixelTypeAsString(pixelType));

    // IOComponentType is UCHAR, CHAR, USHORT, SHORT, UINT, INT, ULONG, LONG, FLOAT, DOUBLE
    itk::ImageIOBase::IOComponentType componentType = imageIO->GetComponentType();
    LINFO("ComponentType: " << imageIO->GetComponentTypeAsString(componentType));

    // NumberOfComponents is usually one, but for non-scalar pixel types, it can be anything
    unsigned int NumberOfComponents = imageIO->GetNumberOfComponents();
    LINFO("Number of Components: " << NumberOfComponents);
    if(NumberOfComponents != 1) {
        throw tgt::UnsupportedFormatException("Unsupported number of components!");
    }

    //-------Info we don't need here:---------------
    //unsigned dims[32];   // almost always no more than 4 dims, but ...
    //unsigned origin[32];
    double spacing[32];
    //std::vector<double> directions[32];
    for(unsigned i = 0; i < NumberOfDimensions && i < 32; i++)
    {
        //dims[i] = imageIO->GetDimensions(i);
        //origin[i] = imageIO->GetOrigin(i);
        spacing[i] = imageIO->GetSpacing(i);
        //directions[i] = imageIO->GetDirection(i);
    }

    Volume* dataset;
    switch(pixelType) {
        case itk::ImageIOBase::SCALAR:
            switch(componentType) {
                case itk::ImageIOBase::UCHAR:
                    dataset = readScalarVolume<uint8_t>(fileName);
                    break;
                case itk::ImageIOBase::CHAR:
                    dataset = readScalarVolume<int8_t>(fileName);
                    break;
                case itk::ImageIOBase::USHORT:
                    dataset = readScalarVolume<uint16_t>(fileName);
                    break;
                case itk::ImageIOBase::SHORT:
                    dataset = readScalarVolume<int16_t>(fileName);
                    break;
                case itk::ImageIOBase::UINT:
                    dataset = readScalarVolume<uint32_t>(fileName);
                    break;
                case itk::ImageIOBase::INT:
                    dataset = readScalarVolume<int32_t>(fileName);
                    break;
#ifndef WIN32
                case itk::ImageIOBase::ULONG:
                    dataset = readScalarVolume<uint64_t>(fileName);
                    break;
                case itk::ImageIOBase::LONG:
                    dataset = readScalarVolume<int64_t>(fileName);
                    break;
#endif
                case itk::ImageIOBase::FLOAT:
                    dataset = readScalarVolume<float>(fileName);
                    break;
                case itk::ImageIOBase::DOUBLE:
                    dataset = readScalarVolume<double>(fileName);
                    break;
                default:
                    throw tgt::UnsupportedFormatException("Unsupported component type!");
            }
            break;
        case itk::ImageIOBase::RGB:
        case itk::ImageIOBase::RGBA:
        case itk::ImageIOBase::VECTOR:
        case itk::ImageIOBase::COVARIANTVECTOR:
        case itk::ImageIOBase::POINT:
        default:
            throw tgt::UnsupportedFormatException("Unsupported pixel type!");
        return 0;
    }

    VolumeCollection* volumeCollection = new VolumeCollection();
    dataset->setOrigin(fileName);
    dataset->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping(dataset->getBaseType()));
    volumeCollection->add(dataset);

    return volumeCollection;
}

VolumeReader* ITKVolumeReader::create(ProgressBar* /*progress*/) const {
    return new ITKVolumeReader(/*progress*/);
}

} // namespace voreen
