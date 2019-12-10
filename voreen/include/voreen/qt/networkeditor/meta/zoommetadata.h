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

#ifndef VRN_ZOOMMETADATA_H
#define VRN_ZOOMMETADATA_H

#include "voreen/core/io/serialization/serialization.h"

#include "voreen/qt/voreenqtapi.h"

namespace voreen {

/**
 * The @c ZoomMetaData class stores some zoom properties like zoom factor or focus position.
 *
 * @see MetaDataBase
 */

class VRN_QT_API ZoomMetaData : public MetaDataBase {
public:
    ZoomMetaData(const tgt::Matrix3d transform = tgt::Matrix3d());
    virtual ~ZoomMetaData();

    virtual std::string getClassName() const { return "ZoomMetaData"; }
    virtual MetaDataBase* create() const;
    virtual MetaDataBase* clone() const;
    virtual std::string toString() const;
    virtual std::string toString(const std::string& component) const;

    virtual void serialize(Serializer& s) const; ///< @see Serializable::serialize
    virtual void deserialize(Deserializer& s);  ///< @see Serializable::deserialize
    void setTransform(tgt::mat3 transform);
    tgt::Matrix3d getTransform() const;

private:
    tgt::Matrix3d transform_;
};


} // namespace

#endif // VRN_ZOOMMETADATA_H
