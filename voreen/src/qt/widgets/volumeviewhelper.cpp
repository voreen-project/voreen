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

#include "voreen/qt/widgets/volumeviewhelper.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumepreview.h"

#include "tgt/filesystem.h"

#include <QDir>

namespace voreen {

std::string VolumeViewHelper::getVolumeType(const VolumeBase* volume) {
    return volume->getFormat();
}

QPixmap VolumeViewHelper::generatePreview(const VolumeBase* volume, int height) {
    return generateBorderedPreview(volume, height, 0);
}

QPixmap VolumeViewHelper::generateBorderedPreview(const VolumeBase* handle, int height, int border) {

    VolumePreview* prev = handle->getDerivedData<VolumePreview>();
    if (prev) {

        int internHeight = prev->getHeight();

        QImage origImg = QImage(internHeight, internHeight, QImage::Format_ARGB32);
        for (int y=0; y<internHeight; y++) {
            for (int x=0; x<internHeight; x++) {
                int previewIndex = y * internHeight + x;
                int greyVal = prev->getData()[previewIndex];
                origImg.setPixel(x, y, qRgb(greyVal, greyVal, greyVal));
            }
        }
        QImage* preview = new QImage(origImg.scaled(height, height, Qt::IgnoreAspectRatio, Qt::SmoothTransformation));

        // draw border
        for (int y = 0; y < preview->height(); y++){
            for (int x = 0; x < border; x++)
                preview->setPixel(preview->width()-x-1, y, qRgb(255, 255, 255));
        }
        for (int x = 0; x < preview->width(); x++){
            for (int y = 0; y < border; y++)
                   preview->setPixel(x, preview->height()-y-1, qRgb(255, 255, 255));
        }

        QPixmap pixmap = QPixmap::fromImage(*preview);
        delete preview;

        return pixmap;
    }
    else {
        LWARNINGC("voreenqt.VolumeViewHelper", "VolumePreview not available");
        return QPixmap(height, height);
    }
}

std::string VolumeViewHelper::volumeInfoString(const VolumeBase* handle) {
    std::string outString;
    if (handle && handle->getRepresentation<VolumeRAM>()) {

        std::string filename = getStrippedVolumeName(handle);
        if(filename.length() > 13) {
            QString ext = QString::fromStdString(filename);
            filename.resize(13);
            filename+="...";
            filename+=ext.section('.',-1).toStdString();
        }
        std::string spacing = getVolumeSpacing(handle);
        std::string dimension = getVolumeDimension(handle);
        std::string type = getVolumeType(handle);
        outString+="<style type='text/css'> table {margin-left:0px; margin-top:1px; font-size: 9px;}</style><table><tr><td>"+filename+"</td><td></td></tr><tr><td>type </td><td>" \
            +type+"</td></tr><tr><td>dimension </td><td>" \
            +dimension+"</td></tr><tr><td>spacing </td><td>" \
            +spacing+"</td></tr><tr><td>memSize </td><td>" \
    +getVolumeMemorySize(handle)+" bytes</td></tr> </table>";

    }
    else {
        outString = "<strong>Volume:</strong> no Volume";
    }
    return outString;
}

std::string VolumeViewHelper::getVolumeName(const VolumeBase* handle){
    std::string volumeName;
    if (handle) {
        //const Volume* vh = dynamic_cast<const Volume*>(handle);
        //if(vh)
        //    volumeName = vh->getOrigin().getPath();
        volumeName = handle->getOrigin().getPath();
    }
    else {
        volumeName = "no Volume";
    }
    return volumeName;
}

std::string VolumeViewHelper::getStrippedVolumeName(const VolumeBase* handle) {
    if (!handle)
        return "no volume";
    else {
        std::string volumePath = getVolumeName(handle);

        //  remove trailing slash (dicom)
        std::string::size_type separator = volumePath.find_last_of("/\\");
        if (separator == volumePath.size()-1)
            volumePath = volumePath.substr(0, volumePath.size()-1);

        return tgt::FileSystem::fileName(volumePath);
    }
}

std::string VolumeViewHelper::getVolumePath(const VolumeBase* handle) {
    std::string volumeName;
    if (handle) {
        //const Volume* vh = dynamic_cast<const Volume*>(handle);
        //if(vh)
        //    volumeName = tgt::FileSystem::dirName(handle->getOrigin().getPath());
        volumeName = handle->getOrigin().getPath();
    }
    else {
        volumeName = "no Volume";
    }
    return volumeName;
}

std::string VolumeViewHelper::getVolumeDimension(const VolumeBase* volume) {
    std::stringstream out;

    out << volume->getDimensions()[0] << " x " << volume->getDimensions()[1] << " x " << volume->getDimensions()[2];

    return out.str();
}

std::string VolumeViewHelper::getVolumeSpacing(const VolumeBase* volume) {
    std::stringstream out;
    out << volume->getSpacing()[0] << " x " << volume->getSpacing()[1] << " x " << volume->getSpacing()[2] << " mm";
    return out.str();
}

std::string VolumeViewHelper::getVolumeMemorySize(const VolumeBase* volume) {
    return formatMemorySize(getVolumeMemorySizeByte(volume));
}

uint64_t VolumeViewHelper::getVolumeMemorySizeByte(const VolumeBase* volume) {
    return static_cast<uint64_t>(volume->getNumVoxels()) * volume->getBytesPerVoxel();

}

} // namespace
