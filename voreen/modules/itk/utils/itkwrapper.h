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

#ifndef VRN_ITKWRAPPER_H
#define VRN_ITKWRAPPER_H

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include <itkImage.h>
#include <itkImportImageFilter.h>
#include "itkImageRegionConstIterator.h"
#include "itkRGBPixel.h"

namespace voreen {

void transferTransformation(const VolumeBase* input, Volume* output);
void transferRWM(const VolumeBase* input, Volume* output);

template<class T>
typename itk::Image<T, 3>::Pointer voreenToITK(const VolumeBase* handle) {
    const VolumeAtomic<T>* vol = dynamic_cast<const VolumeAtomic<T>*>(handle->getRepresentation<VolumeRAM>());
    typedef typename itk::Image<T, 3> ImageType;
    typedef typename ImageType::PixelType PixelType;
    typedef itk::ImportImageFilter<T, 3> ImportFilterType;
    typename ImportFilterType::Pointer importer = ImportFilterType::New();

    typename itk::Image<T, 3>::SizeType size;
    size[0] = vol->getDimensions().x;
    size[1] = vol->getDimensions().y;
    size[2] = vol->getDimensions().z;

    typename itk::Image<T, 3>::IndexType start;
    start.Fill(0);

    typename itk::Image<T, 3>::RegionType region;
    region.SetSize(size);
    region.SetIndex(start);
    importer->SetRegion(region);

    double spacing[3];
    spacing[0] = handle->getSpacing().x;
    spacing[1] = handle->getSpacing().y;
    spacing[2] = handle->getSpacing().z;
    importer->SetSpacing(spacing);

    double origin[3];
    origin[0] = handle->getOffset().x;
    origin[1] = handle->getOffset().y;
    origin[2] = handle->getOffset().z;
    importer->SetOrigin(origin);

    const bool importFilterWillDeleteTheInputBuffer = false;
    PixelType* pixelData = static_cast<PixelType*>(const_cast<void*>(vol->getData())); //ugly...but we do not modify the input data
    const unsigned int totalNumberOfPixels = vol->getNumVoxels();
    importer->SetImportPointer(pixelData, totalNumberOfPixels, importFilterWillDeleteTheInputBuffer);
    importer->Update();

    return importer->GetOutput();
}

// 2d Vector
template<class T>
typename itk::Image<itk::CovariantVector<T,2>, 3>::Pointer voreenVec2ToITKVec2(const VolumeBase* handle) {
    const VolumeAtomic<tgt::Vector2<T> >* vol = dynamic_cast<const VolumeAtomic<tgt::Vector2<T> >*>(handle->getRepresentation<VolumeRAM>());
    typedef typename itk::Image<itk::CovariantVector<T,2>, 3> ImageType;
    typedef typename ImageType::PixelType PixelType;
    typedef itk::ImportImageFilter<itk::CovariantVector<T,2>,3> ImportFilterType;
    typename ImportFilterType::Pointer importer = ImportFilterType::New();

    typename itk::Image<itk::CovariantVector<T,2>, 3>::SizeType size;
    size[0] = vol->getDimensions().x;
    size[1] = vol->getDimensions().y;
    size[2] = vol->getDimensions().z;

    typename itk::Image<itk::CovariantVector<T,2>, 3>::IndexType start;
    start.Fill(0);

    typename itk::Image<itk::CovariantVector<T,2>, 3>::RegionType region;
    region.SetSize(size);
    region.SetIndex(start);
    importer->SetRegion(region);

    double spacing[3];
    spacing[0] = handle->getSpacing().x;
    spacing[1] = handle->getSpacing().y;
    spacing[2] = handle->getSpacing().z;
    importer->SetSpacing(spacing);

    double origin[3];
    origin[0] = handle->getOffset().x;
    origin[1] = handle->getOffset().y;
    origin[2] = handle->getOffset().z;
    importer->SetOrigin(origin);

    const bool importFilterWillDeleteTheInputBuffer = false;
    PixelType* pixelData = static_cast<PixelType*>(const_cast<void*>(vol->getData())); //ugly...but we do not modify the input data
    const unsigned int totalNumberOfPixels = vol->getNumVoxels();
    importer->SetImportPointer(pixelData, totalNumberOfPixels, importFilterWillDeleteTheInputBuffer);
    importer->Update();

    return importer->GetOutput();
}

//3d vector
template<class T>
typename itk::Image<itk::CovariantVector<T,3>, 3>::Pointer voreenVec3ToITKVec3(const VolumeBase* handle) {
    const VolumeAtomic<tgt::Vector3<T> >* vol = dynamic_cast<const VolumeAtomic<tgt::Vector3<T> >*>(handle->getRepresentation<VolumeRAM>());
    typedef typename itk::Image<itk::CovariantVector<T,3>, 3> ImageType;
    typedef typename ImageType::PixelType PixelType;
    typedef itk::ImportImageFilter<itk::CovariantVector<T,3>,3> ImportFilterType;
    typename ImportFilterType::Pointer importer = ImportFilterType::New();

    typename itk::Image<itk::CovariantVector<T,3>, 3>::SizeType size;
    size[0] = vol->getDimensions().x;
    size[1] = vol->getDimensions().y;
    size[2] = vol->getDimensions().z;

    typename itk::Image<itk::CovariantVector<T,3>, 3>::IndexType start;
    start.Fill(0);

    typename itk::Image<itk::CovariantVector<T,3>, 3>::RegionType region;
    region.SetSize(size);
    region.SetIndex(start);
    importer->SetRegion(region);

    double spacing[3];
    spacing[0] = handle->getSpacing().x;
    spacing[1] = handle->getSpacing().y;
    spacing[2] = handle->getSpacing().z;
    importer->SetSpacing(spacing);

    double origin[3];
    origin[0] = handle->getOffset().x;
    origin[1] = handle->getOffset().y;
    origin[2] = handle->getOffset().z;
    importer->SetOrigin(origin);

    const bool importFilterWillDeleteTheInputBuffer = false;
    PixelType* pixelData = static_cast<PixelType*>(const_cast<void*>(vol->getData())); //ugly...but we do not modify the input data
    const unsigned int totalNumberOfPixels = vol->getNumVoxels();
    importer->SetImportPointer(pixelData, totalNumberOfPixels, importFilterWillDeleteTheInputBuffer);
    importer->Update();

    return importer->GetOutput();
}

//4d vector
template<class T>
typename itk::Image<itk::CovariantVector<T,4>, 3>::Pointer voreenVec4ToITKVec4(const VolumeBase* handle) {
    const VolumeAtomic<tgt::Vector4<T> >* vol = dynamic_cast<const VolumeAtomic<tgt::Vector4<T> >*>(handle->getRepresentation<VolumeRAM>());
    typedef typename itk::Image<itk::CovariantVector<T,4>, 3> ImageType;
    typedef typename ImageType::PixelType PixelType;
    typedef itk::ImportImageFilter<itk::CovariantVector<T,4>,3> ImportFilterType;
    typename ImportFilterType::Pointer importer = ImportFilterType::New();

    typename itk::Image<itk::CovariantVector<T,4>, 3>::SizeType size;
    size[0] = vol->getDimensions().x;
    size[1] = vol->getDimensions().y;
    size[2] = vol->getDimensions().z;

    typename itk::Image<itk::CovariantVector<T,4>, 3>::IndexType start;
    start.Fill(0);

    typename itk::Image<itk::CovariantVector<T,4>, 3>::RegionType region;
    region.SetSize(size);
    region.SetIndex(start);
    importer->SetRegion(region);

    double spacing[3];
    spacing[0] = handle->getSpacing().x;
    spacing[1] = handle->getSpacing().y;
    spacing[2] = handle->getSpacing().z;
    importer->SetSpacing(spacing);

    double origin[3];
    origin[0] = handle->getOffset().x;
    origin[1] = handle->getOffset().y;
    origin[2] = handle->getOffset().z;
    importer->SetOrigin(origin);

    const bool importFilterWillDeleteTheInputBuffer = false;
    PixelType* pixelData = static_cast<PixelType*>(const_cast<void*>(vol->getData())); //ugly...but we do not modify the input data
    const unsigned int totalNumberOfPixels = vol->getNumVoxels();
    importer->SetImportPointer(pixelData, totalNumberOfPixels, importFilterWillDeleteTheInputBuffer);
    importer->Update();

    return importer->GetOutput();
}

template<class T>
//VolumeAtomic<T>* ITKToVoreenCopy(const typename itk::Image<T, 3>::Pointer vol) {
Volume* ITKToVoreenCopy(const typename itk::Image<T, 3>* vol) {
    typedef typename itk::Image<T, 3> ImageType;
    //typedef typename ImageType::PixelType PixelType;

    typename ImageType::RegionType region;
    region = vol->GetBufferedRegion();

    typename ImageType::SizeType size = region.GetSize();
//    typename ImageType::IndexType start = region.GetIndex();
    tgt::ivec3 dim(size[0], size[1], size[2]);

    VolumeAtomic<T>* out = new VolumeAtomic<T>(dim);

    typename ImageType::SpacingType sp = vol->GetSpacing();
    tgt::vec3 spacing(sp[0], sp[1], sp[2]);

    typename ImageType::PointType o = vol->GetOrigin();
    tgt::vec3 offset(o[0], o[1], o[2]);

    typename itk::ImageRegionConstIterator< ImageType > it(vol, region);
    it.GoToBegin();
    T* data = (T*)out->getData();

    while( ! it.IsAtEnd() )
    {
        *data = it.Get();
        ++it;
        ++data;
    }

    return new Volume(out, spacing, offset);
}

//2d vector
template<class T>
Volume* ITKVec2ToVoreenVec2Copy(const typename itk::Image<itk::CovariantVector<T,2>, 3>* vol) {
    typedef typename itk::Image<itk::CovariantVector<T,2>, 3> ImageType;
    //typedef typename ImageType::PixelType PixelType;

    typename ImageType::RegionType region;
    region = vol->GetBufferedRegion();

    typename ImageType::SizeType size = region.GetSize();
//    typename ImageType::IndexType start = region.GetIndex();
    tgt::ivec3 dim(size[0], size[1], size[2]);

    VolumeAtomic<tgt::Vector2<T> >* out = new VolumeAtomic<tgt::Vector2<T> >(dim);

    typename ImageType::SpacingType sp = vol->GetSpacing();
    tgt::vec3 spacing(sp[0], sp[1], sp[2]);

    typename ImageType::PointType o = vol->GetOrigin();
    tgt::vec3 offset(o[0], o[1], o[2]);

    typename itk::ImageRegionConstIterator< ImageType > it(vol, region);
    it.GoToBegin();
    tgt::col3* data = (tgt::col3*)out->getData();

    while( ! it.IsAtEnd() )
    {
        (*data)[0] = it.Get()[0];
        (*data)[1] = it.Get()[1];
        ++it;
        ++data;
    }

    return new Volume(out, spacing, offset);
}

//3d vector
template<class T>
Volume* ITKVec3ToVoreenVec3Copy(const typename itk::Image<itk::CovariantVector<T,3>, 3>* vol) {
    typedef typename itk::Image<itk::CovariantVector<T,3>, 3> ImageType;
    //typedef typename ImageType::PixelType PixelType;

    typename ImageType::RegionType region;
    region = vol->GetBufferedRegion();

    typename ImageType::SizeType size = region.GetSize();
//    typename ImageType::IndexType start = region.GetIndex();
    tgt::ivec3 dim(size[0], size[1], size[2]);

    VolumeAtomic<tgt::Vector3<T> >* out = new VolumeAtomic<tgt::Vector3<T> >(dim);

    typename ImageType::SpacingType sp = vol->GetSpacing();
    tgt::vec3 spacing(sp[0], sp[1], sp[2]);

    typename ImageType::PointType o = vol->GetOrigin();
    tgt::vec3 offset(o[0], o[1], o[2]);

    typename itk::ImageRegionConstIterator< ImageType > it(vol, region);
    it.GoToBegin();
    tgt::col3* data = (tgt::col3*)out->getData();

    while( ! it.IsAtEnd() )
    {
        (*data)[0] = it.Get()[0];
        (*data)[1] = it.Get()[1];
        (*data)[2] = it.Get()[2];
        ++it;
        ++data;
    }

    return new Volume(out, spacing, offset);
}

//3d vector
template<class T>
Volume* ITKRGBToVoreenVec3Copy(const typename itk::Image<itk::RGBPixel<T>, 3>* vol) {
    typedef typename itk::Image<itk::RGBPixel<T>, 3> ImageType;
    //typedef typename ImageType::PixelType PixelType;

    typename ImageType::RegionType region;
    region = vol->GetBufferedRegion();

    typename ImageType::SizeType size = region.GetSize();
//    typename ImageType::IndexType start = region.GetIndex();
    tgt::ivec3 dim(size[0], size[1], size[2]);

    VolumeAtomic<tgt::Vector3<T> >* out = new VolumeAtomic<tgt::Vector3<T> >(dim);

    typename ImageType::SpacingType sp = vol->GetSpacing();
    tgt::vec3 spacing(sp[0], sp[1], sp[2]);

    typename ImageType::PointType o = vol->GetOrigin();
    tgt::vec3 offset(o[0], o[1], o[2]);

    typename itk::ImageRegionConstIterator< ImageType > it(vol, region);
    it.GoToBegin();
    tgt::col3* data = (tgt::col3*)out->getData();

    while( ! it.IsAtEnd() )
    {
        (*data)[0] = it.Get()[0];
        (*data)[1] = it.Get()[1];
        (*data)[2] = it.Get()[2];
        ++it;
        ++data;
    }

    return new Volume(out, spacing, offset);
}

//4d vector
template<class T>
Volume* ITKVec4ToVoreenVec4Copy(const typename itk::Image<itk::CovariantVector<T,4>, 3>* vol) {
    typedef typename itk::Image<itk::CovariantVector<T,4>, 3> ImageType;
    //typedef typename ImageType::PixelType PixelType;

    typename ImageType::RegionType region;
    region = vol->GetBufferedRegion();

    typename ImageType::SizeType size = region.GetSize();
//    typename ImageType::IndexType start = region.GetIndex();
    tgt::ivec3 dim(size[0], size[1], size[2]);

    VolumeAtomic<tgt::Vector4<T> >* out = new VolumeAtomic<tgt::Vector4<T> >(dim);

    typename ImageType::SpacingType sp = vol->GetSpacing();
    tgt::vec3 spacing(sp[0], sp[1], sp[2]);

    typename ImageType::PointType o = vol->GetOrigin();
    tgt::vec3 offset(o[0], o[1], o[2]);

    typename itk::ImageRegionConstIterator< ImageType > it(vol, region);
    it.GoToBegin();
    tgt::col3* data = (tgt::col3*)out->getData();

    while( ! it.IsAtEnd() )
    {
        (*data)[0] = it.Get()[0];
        (*data)[1] = it.Get()[1];
        (*data)[2] = it.Get()[2];
        (*data)[3] = it.Get()[3];
        ++it;
        ++data;
    }

    return new Volume(out, spacing, offset);
}

}   //namespace

#endif
