/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#ifndef VRN_TRACKBALL_H
#define VRN_TRACKBALL_H

#include <cmath>
#include <vector>
#include "voreen/core/voreencoreapi.h"
#include "voreen/core/voreenobject.h"
#include "tgt/camera.h"

namespace voreen {

class CameraProperty;

/**
 *  This class implements a Trackball which can be used to freely rotate an object
 *  around a given center.
 */
class VRN_CORE_API VoreenTrackball : public VoreenSerializableObject {

public:
    VoreenTrackball(CameraProperty* camera);
    virtual ~VoreenTrackball();
    virtual VoreenSerializableObject* create() const;

    virtual std::string getClassName() const  { return "VoreenTrackball"; };

    /// Resets the trackball to the initial configuration of the canvas' camera.
    //void reset();

    /// Returns the camera this trackball operates on.
    CameraProperty* getCamera() const;

    /// Assigns the camera to be modified by the trackball.
    void setCamera(CameraProperty* camera);

    /// Rotate the trackball according to Quaternion quat.
    /// @param quat Quaternion represention rotation.
    void rotate(tgt::Quaternion<float> quat);

    /// Rotate the trackball by angle phi around axis axis.
    /// @param axis axis in camera coordinates.
    virtual void rotate(tgt::vec3 axis, float phi);

    /// Rotate the trackball according to new mouse position.
    /// @param mouse coodinates of mouse, scaled to [-1, 1]x[-1, 1]
    void rotate(const tgt::vec2& newMouse, const tgt::vec2& lastMouse);

    /// Move the trackball along axis.
    /// @param length Distance relative to frustum dimensions at trackball center.
    ///               Trackball center will get out of sight when >> 1.
    /// @param axis Axis in camera coordinates along which to move.
    virtual void move(float length, tgt::vec3 axis);

    /// Move the trackball according to new mouse position.
    /// @param mouse coodinates of mouse, scaled to [-1, 1]x[-1, 1]
    void move(const tgt::vec2& newMouse, const tgt::vec2& lastMouse);

    /// Zoom in by factor factor.
    virtual void zoom(float factor);

    /// Zoom in according to new mouse position.
    /// @param mouse coodinates of mouse, scaled to [-1, 1]x[-1, 1]
    void zoom(const tgt::vec2& newMouse, const tgt::vec2& lastMouse, const tgt::vec2& mouseZoomInDirection);

    /// set an absolute Distance from Focus
    void zoomAbsolute(float focallength);

    /// sets the ortho zoom factor according to the focal length
    void adjustOrthogonalZoomToFocalLength();

    /// last rotation applied to trackball
    tgt::quat getLastOrientationChange() const { return lastOrientationChange_; }

    /// If the center of the trackball should be moved relative to objects.
    /// This has large influence on trackball behavior. If you look at one object located at
    /// world coordinates origin, you might want the center not to move. This is used in voreen.
    /// If you use trackball to look at mutiple objects, located at varying position, you might
    /// want the center to move.
    /// See trackball demo in samples folder.
    void setMoveCenter(bool b = true) { moveCenter_ = b; };
    bool getMoveCenter() const { return moveCenter_; };

    /// Set trackballs center. Most people do not need this method. Use with caution!
    void setCenter(tgt::vec3 center) { center_ = center; };
    tgt::vec3 getCenter() const { return center_; };

    /// Move the camera in a certain direction
    void moveCamera(const tgt::vec3& motionvector);

    /// set the radius of the trackball to a given value
    void setSize(float size) { size_ = size; };
    float getSize() const { return size_; };

    /// Returns the orthogonal distance between camera position and center of the trackball
    /// with respect to the camera's look vector.
    float getCenterDistance();

    virtual void serialize(Serializer& s) const;

    virtual void deserialize(Deserializer& s);

protected:

    CameraProperty* cameraProperty_;

    /// this holds the center around which the camera will be rotated
    tgt::vec3 center_;

    /// Wheather to change the center in move method.
    /// center_ will alway be (0,0,0) if false.
    bool moveCenter_;

    /// the trackball size effects the trackball sensitivity
    float size_;

    /// store the initial position, focus point and up vector of camera used to be able to
    /// reset camera to those values
    /*void saveCameraParameters();
    vec3 cameraPosition_;
    vec3 cameraFocus_;
    vec3 cameraUpVector_; */

    /// last rotation applied to trackball
    tgt::quat lastOrientationChange_;

    /// scale screen-coodinates of mouse to intervall [-1, 1]x[-1, 1]
    tgt::vec2 scaleMouse(const tgt::ivec2& sc) const;

    /// projection math
    tgt::vec3 projectToSphere(const tgt::vec2& xy) const;

    /// transform vectors given in camera coordinates to vectors in world coordinates
    tgt::vec3 coordTransform(const tgt::vec3& axis) const;

    float getRotationAngle(float acuteness) const;
    float getMovementLength(float acuteness) const;
    float getZoomFactor(float acuteness, bool zoomIn) const;
    float getRollAngle(float acuteness, bool left) const;

};

} // namespace voreen

#endif // VRN_TRACKBALL_H
