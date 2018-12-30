/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "dicominfo.h"

DicomInfo::DicomInfo()
{
}

void DicomInfo::setDx(int dx){
    dx_ = dx;
}

void DicomInfo::setDy(int dy){
    dy_ = dy;
}

void DicomInfo::setDz(int dz){
    dz_ = dz;
}

void DicomInfo::setNumberOfFrames(int frames) {
    numberOfFrames_ = frames;
}

void DicomInfo::setXSpacing(double x){
    xSpacing_ = x;
}

void DicomInfo::setYSpacing(double y){
    ySpacing_ = y;
}

void DicomInfo::setZSpacing(double z){
    zSpacing_ = z;
}

void DicomInfo::setXOrientationPatient(tgt::dvec3 orientation){
    xOrientationPatient_ = orientation;
}

void DicomInfo::setYOrientationPatient(tgt::dvec3 orientation){
    yOrientationPatient_ = orientation;
}

void DicomInfo::setSliceNormal(tgt::dvec3 normal){
    sliceNormal_ = normal;
}

void DicomInfo::setOffset(tgt::dvec3 offset){
    offset_ = offset;
}

void DicomInfo::setBitsStored(int bits){
    bitsStored_ = bits;
}

void DicomInfo::setSamplesPerPixel(int spp){
    samplesPerPixel_ = spp;
}

void DicomInfo::setBytesPerVoxel(int bytes){
    bytesPerVoxel_ = bytes;
}

void DicomInfo::setPixelRepresentation(unsigned short representation) {
    pixelRepresentation_ = representation;
}

void DicomInfo::setIntercept(float intercept){
    intercept_ = intercept;
}

void DicomInfo::setSlope(float slope){
    slope_ = slope;
}

void DicomInfo::setRescaleType(std::string type){
    rescaleType_ = type;
}

void DicomInfo::setModality(std::string modality){
    modality_ = modality;
}

void DicomInfo::setSeriesInstanceUID(std::string uid){
    seriesInstanceUID_ = uid;
}

void DicomInfo::setStudyInstanceUID(std::string uid){
    studyInstanceUID_ = uid;
}

void DicomInfo::setSeriesDescription(std::string description){
    seriesDescription_ = description;
}

void DicomInfo::setStudyDescription(std::string description){
    studyDescription_ = description;
}

void DicomInfo::setPatientName(std::string name){
    patientName_ = name;
}

void DicomInfo::setPatientId(std::string id){
    patientId_ = id;
}

void DicomInfo::setRwmDiffers(bool d) {
    rwmDiffers_ = d;
}

void DicomInfo::setFormat(std::string format) {
    format_ = format;
}

void DicomInfo::setBaseType(std::string type) {
    baseType_ = type;
}

int DicomInfo::getDx() const {
    return dx_;
}

int DicomInfo::getDy() const {
    return dy_;
}


int DicomInfo::getDz() const {
    return dz_;
}

int DicomInfo::getNumberOfFrames() const {
    return numberOfFrames_;
}

double DicomInfo::getXSpacing() const {
    return xSpacing_;
}

double DicomInfo::getYSpacing() const {
    return ySpacing_;
}

double DicomInfo::getZSpacing() const {
    return zSpacing_;
}

tgt::dvec3 DicomInfo::getXOrientationPatient() const {
    return xOrientationPatient_;
}

tgt::dvec3 DicomInfo::getYOrientationPatient() const {
    return yOrientationPatient_;
}

tgt::dvec3 DicomInfo::getSliceNormal() const {
    return sliceNormal_;
}

tgt::dvec3 DicomInfo::getOffset() const {
    return offset_;
}

int DicomInfo::getBitsStored() const {
    return bitsStored_;
}

int DicomInfo::getSamplesPerPixel() const {
    return samplesPerPixel_;
}

unsigned short DicomInfo::getPixelRepresentation() const {
    return pixelRepresentation_;
}

int DicomInfo::getBytesPerVoxel() const {
    return bytesPerVoxel_;
}

float DicomInfo::getIntercept() const {
    return intercept_;
}

float DicomInfo::getSlope() const {
    return slope_;
}

std::string DicomInfo::getRescaleType() const {
    return rescaleType_;
}

std::string DicomInfo::getModality() const {
    return modality_;
}

std::string DicomInfo::getSeriesInstanceUID() const {
    return seriesInstanceUID_;
}

std::string DicomInfo::getStudyInstanceUID() const {
    return studyInstanceUID_;
}

std::string DicomInfo::getSeriesDescription() const {
    return seriesDescription_;
}

std::string DicomInfo::getStudyDescription() const {
    return studyDescription_;
}

std::string DicomInfo::getPatientName() const {
    return patientName_;
}

std::string DicomInfo::getPatientId() const {
    return patientId_;
}

bool DicomInfo::rwmDiffers() const {
    return rwmDiffers_;
}

std::string DicomInfo::getFormat() const {
    return format_;
}

std::string DicomInfo::getBaseType() const {
    return baseType_;
}
