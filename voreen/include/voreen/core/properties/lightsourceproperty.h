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

#ifndef VRN_LIGHTSOURCEPROPERTY_H
#define VRN_LIGHTSOURCEPROPERTY_H

#include "voreen/core/properties/vectorproperty.h"
#include "tgt/camera.h"

namespace voreen {

class CameraProperty;

class VRN_CORE_API LightSourceProperty : public FloatVec4Property {

public:

    LightSourceProperty();

    LightSourceProperty(const std::string& id, const std::string& guiText, const tgt::vec4& value,
                  const tgt::vec4& minimum = tgt::vec4(-10000000.0f), const tgt::vec4& maximum = tgt::vec4(10000000.0f),
                  int invalidationLevel=Processor::INVALID_RESULT, Property::LevelOfDetail lod = Property::LOD_DEFAULT);

    virtual Property* create() const {
        return new LightSourceProperty();
    }

    virtual std::string getClassName() const       { return "LightSourceProperty"; }
    virtual std::string getTypeDescription() const { return "LightSource"; }

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    void setLightPos(const tgt::vec4& lightPos);
    tgt::vec4 getLightPos() const;

    void setFollowCam(bool);
    bool getFollowCam() const;

    void setMaxDist(float d);
    float getMaxDist() const;

    bool hasCamera() const;
    CameraProperty* getCamera();
    void setCamera(CameraProperty* cam);

    void set(const tgt::vec4& value);

protected:

    bool followCam_;
    tgt::vec3 curCenter_;
    float maxDist_;
    tgt::vec4 lightPos_;
    CameraProperty* camProp_;

    void cameraUpdate();    // on property change of the associated cameraproperty this is invoked by the onChange mechanism
};

} // namespace voreen

#endif // VRN_LIGHTSOURCEPROPERTY_H
