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

#include "voreen/qt/networkeditor/meta/zoommetadata.h"

namespace voreen {

ZoomMetaData::ZoomMetaData(const tgt::Matrix3d transform)
    : transform_(transform)
{}

ZoomMetaData::~ZoomMetaData() {}

void ZoomMetaData::serialize(Serializer& s) const {
    s.serialize("transform", transform_);
}

void ZoomMetaData::deserialize(Deserializer& s) {
    s.deserialize("transform", transform_);
}

void ZoomMetaData::setTransform(tgt::mat3 transform) {
    transform_ = transform;
}

tgt::Matrix3d ZoomMetaData::getTransform() const {
    return transform_;
}

MetaDataBase* ZoomMetaData::clone() const {
    return new ZoomMetaData(transform_);
}

MetaDataBase* ZoomMetaData::create() const {
    return new ZoomMetaData();
}

std::string ZoomMetaData::toString() const {
    std::stringstream s;
    s << transform_;
    return s.str();
}

std::string ZoomMetaData::toString(const std::string& /*component*/) const {
    return toString();
}

} // namespace
