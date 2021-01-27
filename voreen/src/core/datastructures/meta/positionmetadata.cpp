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

#include "voreen/core/datastructures/meta/positionmetadata.h"
#include "voreen/core/utils/stringutils.h"

namespace voreen {

PositionMetaData::PositionMetaData(const int& x, const int& y)
    : x_(x)
    , y_(y)
{}

PositionMetaData::~PositionMetaData() {}

void PositionMetaData::serialize(Serializer& s) const {
    s.serialize("x", x_);
    s.serialize("y", y_);
}

void PositionMetaData::deserialize(Deserializer& s) {
    s.deserialize("x", x_);
    s.deserialize("y", y_);
}

void PositionMetaData::setX(const int& value) {
    x_ = value;
}

int PositionMetaData::getX() const {
    return x_;
}

void PositionMetaData::setY(const int& value) {
    y_ = value;
}

int PositionMetaData::getY() const {
    return y_;
}

MetaDataBase* PositionMetaData::clone() const {
    return new PositionMetaData(x_, y_);
}

MetaDataBase* PositionMetaData::create() const {
    return new PositionMetaData();
}

std::string PositionMetaData::toString() const {
    return (itos(x_) + ", " + itos(y_));
}

std::string PositionMetaData::toString(const std::string& component) const {
    if (component == "x")
        return itos(x_);
    else if (component == "y")
        return itos(y_);
    else
        return toString();
}

} // namespace
