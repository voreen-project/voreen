/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_CAMERAPROPERTY_H
#define VRN_CAMERAPROPERTY_H

#include "voreen/core/properties/templateproperty.h"
#include "voreen/core/properties/condition.h"
#include "voreen/core/interaction/voreentrackball.h"
#include "tgt/camera.h"

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API TemplateProperty<tgt::Camera>;
#endif

class VoreenTrackball;

/**
 * Property encapsulating a tgt::Camera object.
 */
class VRN_CORE_API CameraProperty : public TemplateProperty<tgt::Camera> {

public:

    enum TrackballCenter {
        SCENE,
        WORLD,
        CAMSHIFT
    };

    /**
     * Constructor.
     *
     * @param scenedAdjuster When set to true, the cameras property widget will include
     *        options to control the camera's behaviour if a scene size change occurs. This
     *        should only be enabled for processors with volume ports or other means to actually
     *        determine the scene size.
     * @param adjustProjectionToViewport when set to true, the camera's projection matrix
     *        is adjusted automatically to viewport changes. This is especially necessary to
     *        reflect the viewport's aspect ratio.
     */
    CameraProperty(const std::string& id, const std::string& guiText, tgt::Camera const value = tgt::Camera(),
               bool sceneAdjuster = false,
               bool adjustProjectionToViewport = true,
               float maxValue = 500.f,
               int invalidationLevel=Processor::INVALID_RESULT, Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    CameraProperty();
    virtual ~CameraProperty();

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "CameraProperty"; }
    virtual std::string getTypeDescription() const { return "Camera"; }

    /**
     * Assigns the passed camera object to the property.
     */
    void set(const tgt::Camera& camera);
    void setPosition(const tgt::vec3& pos);
    void setFocus(const tgt::vec3& focus);
    void setUpVector(const tgt::vec3& up);
    void setFrustum(const tgt::Frustum& frust);
    void setFarDist(float dist);
    void setNearDist(float dist);

    void setUseOrthoZoomFlag(bool b);
    void setOrthoZoomFactor(float zoomFactor);

    bool setStereoEyeSeparation(float separation);
    bool setStereoFocalLength(float focallength);
    bool setStereoWidth(float width);
    bool setStereoEyeMode(tgt::Camera::StereoEyeMode mode);
    bool setStereoRelativeFocalLength(float stereoRelativeFocalLength);
    bool setUseRealWorldFrustum(bool useRealWorldFrustum);
    //bool stereoShift(tgt::vec3 shift);
    bool setStereoAxisMode(tgt::Camera::StereoAxisMode mode);

    void adaptInteractionToScene(const tgt::Bounds& bounds, float nearDist = 0.f, bool forceUpdate = false);
    void resetCameraFocusToTrackballCenter();

    /// sets the ortho zoom factor according to the focal length
    void adjustOrthogonalZoomToFocalLength();

    void setAdaptOnChange(bool b);
    bool getAdaptOnChange() const;

    void setTrackballCenterBehaviour(TrackballCenter t);
    TrackballCenter getTrackballCenterBehaviour() const;

    /// Set the minimum absolute value that an element of the camera position vector can have
    void setMinValue(float val);

    /// \sa setMinValue
    float getMinValue() const;

    /// Set the maximum absolute value that an element of the camera position vector can have
    void setMaxValue(float val);

    /// \sa setMaxValue
    float getMaxValue() const;

    /**
     * @see Property::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Property::deserialize
     */
    virtual void deserialize(Deserializer& s);

    virtual void look(tgt::ivec2 viewportSize);
    VoreenTrackball& getTrackball();
    const VoreenTrackball& getTrackball() const;

    bool isSceneAdjuster() const;

    const tgt::Bounds& getSceneBounds() const;
    void setSceneBounds(const tgt::Bounds& b);

private:

    bool sceneAdjuster_;
    VoreenTrackball trackball_;
    float minValue_;
    float maxValue_;

    // Stores the current bounds of the scene (if present) to adapt interaction to the scene size
    tgt::Bounds currentSceneBounds_;

    bool adaptOnChange_;
    TrackballCenter centerOption_;
    //bool autoFocus_; //currently not in use
};

} // namespace voreen

#endif // VRN_CAMERAPROPERTY_H
