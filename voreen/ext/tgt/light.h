/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 Visualization and Computer Graphics Group, *
 * Department of Computer Science, University of Muenster, Germany.   *
 * <http://viscg.uni-muenster.de>                                     *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#ifndef TGT_LIGHT_H
#define TGT_LIGHT_H

#include "tgt/vector.h"
#include "tgt/spline.h"
#include "tgt/curve.h"

#include <list>

#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE

namespace tgt {


/// \todo focus point for spotlights
/// \todo what about spotlight direction while moving along space curve?


/**
    Wraps around OpenGL lights. \n
    This class cannot be instantiated directly, you rather create a light source
    object by calling the static function createLight(). The maximum number of
    concurrent lights is limited by OpenGL. This limit is implementation dependend,
    but at least 8 light sources are guaranteed by the OpenGL specification. \n
    The OpenGL light source identifier (GL_LIGHT0 .. GL_LIGHTx) is transparent to the user,
    createLight() always assigns the next free index to the light source it creates. \n
    All setter calls immediately change the OpenGL state - there is no activate() method. \n
    There is also a special feature: You can move light objects along a space curve! \n
    \n
    Additionally this class provides static methods for controlling the global light model settings.\n
    \n
    See <a href="http://www.opengl.org/documentation/specs/man_pages/hardcopy/GL/html/gl/light.html">glLight</a>
    and <a href="http://www.opengl.org/documentation/specs/man_pages/hardcopy/GL/html/gl/lightmodel.html">glLightModel</a>
    documentation for information about OpenGL light sources and global light model features.
 */
class Light {
public:

    //! Creates a light object and returns a pointer to it. \n
    /// If already MAX_LIGHTS light objects are instantiated, this method returns a null-pointer.
    /// You can check wether a light index is available by calling lightLimitReached(). \n
    /// All parameters of the created light are set to OpenGL defaults.
    /// \note{createLight always choses the next free light index starting at GL_LIGHT0}
    static Light *createLight(bool enable = true){

        // find the next free light index
        int index = getNextLightIndex();

        // found free light index -> create new light
        if (index >= GL_LIGHT0){
            Light *light = new Light(index);
            light->setAttributes(vec4(0.f, 0.f, 1.f, 0.f), vec4(0.f, 0.f, 0.f, 1.f), vec4(1.f, 1.f, 1.f, 1.f),
                                 vec4(0.f, 0.f, 0.f, 1.f), vec3(0.f, 0.f, -1.f), 0.f, 180.f, 1.f, 0.f, 0.f);
            light->enable(enable);
            light->setCurve(0);
            usedLightIndices_.push_back(index);
            return light;
        // did not find free light index (light limit reached)
        } else {
            tgtAssert(false,
                "Unable to create new light source. OpenGL light source limit reached.");
            return 0;
        }
    }

    //! Creates a new light source object with the specified parameters and returns a pointer to it.
    /// If light limit is reached, a null-pointer is returned. The unspecified parameters
    /// are set to OpenGL defaults.
    /// \param position Position of the light source. Is transformed to eye-coordinates
    ///                    by the current ModelView-Matrix!
    /// \param ambient,diffuse,specular Light emittance parameters
    /// \param enable specifie wether light source is activated after creation
    static Light *createLight(const vec4& position, const vec4& ambient, const vec4& diffuse, const vec4& specular,
                              bool enable = true){
        Light *light = createLight(enable);
        if (light != 0){
            light->setPosition(position);
            light->setAmbient(ambient);
            light->setDiffuse(diffuse);
            light->setSpecular(specular);
        }

        return light;
    }


    //! Creates a new light source object with the specified parameters and returns a pointer to it.
    /// If light limit is reached, a null-pointer is returned.
    /// \param position Position of the light source. Is transformed to Eye-Coordinates
    ///                          by the current ModelView-Matrix!
    /// \param ambient,diffuse,specular Emittance parameters
    /// \param spotDirection Direction of of the spotlight. This only has an effect,
    ///                          if spotCutoff is != 180. The spotDirection is transformed to
    ///                          eye-coordinates by the current ModelView-Matrix!
    /// \param spotExponent specifies the light distribution of the spotlight. Default value is 0, resulting
    ///                          in uniform light distribution. Only values in the range [0,128] are accepted.
    /// \param spotCutoff specifies the maximum spread angle of a light source. Only values in the
    ///                          range [0,90] and the default value 180 are accepted.
    ///                          All values in range [0,90] create a spotlight. \n
    ///                          See glLight documentation for more information on spotlights.
    /// \param Attenuation Attenuation factors. See glLight documentation for a detailed description of attenuation.
    ///                          Only nonnegative values are accepted.
    /// \param enable specifie wether light source is activated after creation
    static Light *createLight(const vec4& position, const vec4& ambient, const vec4& diffuse, const vec4& specular,
                              const vec3& spotDirection, GLfloat spotExponent, GLfloat spotCutoff,
                              GLfloat constantAttenuation, GLfloat linearAttenuation, GLfloat quadraticAttenuation,
                              bool enable = true){
        Light *light = createLight(position, ambient, diffuse, specular, enable);
        if (light != 0){
            light->setSpotDirection(spotDirection);
            light->setSpotExponent(spotExponent);
            light->setSpotCutoff(spotCutoff);
            light->setAttenuationFactors(constantAttenuation, linearAttenuation, quadraticAttenuation);
        }

        return light;
    }

    // Destructor.
    // Removes the object's lightIndex from the usedLightIndices_ list and disables the light.
    ~Light(){
        usedLightIndices_.remove(lightIndex_);
        enable(false);
    }

    /// Enable or disable light.
    void enable(bool enable = true){
        isEnabled_ = enable;
        if (enable == true){
            glEnable(lightIndex_);
        } else {
            glDisable(lightIndex_);
        }
    }

    /// Check if light is enabled.
    bool isEnabled(){
        return isEnabled_;
    }

    //! Set position of the light source. \n
    /// \note{position is transformed to eye-coordinates by the current ModelView-Matrix!}
    void setPosition(const vec4& position){
        position_ = position;
        glLightfv(lightIndex_, GL_POSITION, position.elem);
    }

    /// Get light source's position as specified by the user (in world-coordinates).
    vec4 getPosition(){
        return position_;
    }

    /// Get light source's position in eye-coordinates (= after ModelView-Transformation)
    vec4 getPositionEye(){
        vec4 pos = vec4();
        glGetLightfv(lightIndex_, GL_POSITION, pos.elem);
        return pos;
    }

    /// Set ambient RGBA intensity of the light source.
    void setAmbient(const vec4& ambient){
        ambient_ = ambient;
        glLightfv(lightIndex_, GL_AMBIENT, ambient.elem);
    }

    /// Get ambient RGBA intensity of the light source.
    vec4 getAmbient(){
        return ambient_;
    }

    /// Set diffuse RGBA intensity of the light source.
    void setDiffuse(const vec4& diffuse){
        diffuse_ = diffuse;
        glLightfv(lightIndex_, GL_DIFFUSE, diffuse.elem);
    }

    /// Get diffuse RGBA intensity of the light source.
    vec4 getDiffuse(){
        return diffuse_;
    }

    /// Set specular RGBA intensity of the light source.
    void setSpecular(const vec4& specular){
        specular_ = specular;
        glLightfv(lightIndex_, GL_SPECULAR, specular.elem);
    }

    /// Get specular RGBA intensity of the light source.
    vec4 getSpecular(){
        return specular_;
    }

    //! Set the light source's spot direction. This only has an effect if the spot cutoff ist < 180. \n
    /// \note{spot direction is transformed to eye-coordinates by the current ModelView-Matrix.}
    void setSpotDirection(const vec3& spotDirection){
        spotDirection_ = spotDirection;
        glLightfv(lightIndex_, GL_SPOT_DIRECTION, spotDirection.elem);
    }


    /// Get light source's spot direction as specified by the user (in world-coordinates).
    vec3 getSpotDirection(){
        return spotDirection_;
    }

    /// Get light source's spot direction in eye-coordinates (= after ModelView-Transformation)
    vec3 getSpotDirectionEye(){
        vec3 pos = vec3();
        glGetLightfv(lightIndex_, GL_SPOT_DIRECTION, pos.elem);
        return pos;
    }

    //! Set light source's spot exponent, which specifies its light distribution.
    /// The default value is 0, resulting in uniform light distribution. \n
    /// Only values in the range [0,128] are accepted.
    void setSpotExponent(GLfloat spotExponent){
        tgtAssert(spotExponent >= 0 && spotExponent <= 128,
            "Invalid Spot Exponent. Must be inside range [0,128].");
        spotExponent_ = spotExponent;
        glLightf(lightIndex_, GL_SPOT_EXPONENT, spotExponent);
    }

    /// Get spot exponent.
    GLfloat getSpotExponent(){
        return spotExponent_;
    }

    //! Set spot cutoff. Specifies the maximum spread angle of the light source. \n
    /// Only values in the range [0,90] and the default value 180 are accepted.
    /// All values in range [0,90] create a spotlight. \n
    ///    See glLight documentation for more information on spotlights.
    void setSpotCutoff(GLfloat spotCutoff){
        tgtAssert((spotCutoff >= 0 && spotCutoff <= 90) || spotCutoff == 180,
            "Invalid Spot Cutoff. Must be inside range [0,90] or special value 180.");
        spotCutoff_ = spotCutoff;
        glLightf(lightIndex_, GL_SPOT_CUTOFF, spotCutoff);
    }

    /// Get spot cutoff.
    GLfloat getSpotCutoff(){
        return spotCutoff_;
    }

    //! Set constant, linear and quadratic attenuation factors. Only nonnegative values are accepted. \n
    /// See glLight documentation for a description of attenuation.
    void setAttenuationFactors(GLfloat constant, GLfloat linear, GLfloat quadratic){
        tgtAssert(constant >= 0 && linear >= 0 && quadratic >= 0,
            "Invalid Attenuation Factor. All factors must be nonnegative.");
        constantAttenuation_ = constant;
        linearAttenuation_ = linear;
        quadraticAttenuation_ = quadratic;
        glLightf(lightIndex_, GL_CONSTANT_ATTENUATION, constant);
        glLightf(lightIndex_, GL_LINEAR_ATTENUATION, linear);
        glLightf(lightIndex_, GL_QUADRATIC_ATTENUATION, quadratic);
    }

    /// Returns the light source's attenuation factors.
    void getAttenuationFactors(GLfloat& constant, GLfloat& linear, GLfloat& quadratic){
        constant = constantAttenuation_;
        linear = linearAttenuation_;
        quadratic = quadraticAttenuation_;
    }

    /// Set all parameters of the light source at once.
    void setAttributes(const vec4& position, const vec4& ambient, const vec4& diffuse, const vec4& specular,
                       const vec3& spotDirection, GLfloat spotExponent, GLfloat spotCutoff,
                       GLfloat constantAttenuation, GLfloat linearAttenuation, GLfloat quadraticAttenuation){
        setPosition(position);
        setAmbient(ambient);
        setDiffuse(diffuse);
        setSpecular(specular);
        setSpotDirection(spotDirection);
        setSpotExponent(spotExponent);
        setSpotCutoff(spotCutoff);
        setAttenuationFactors(constantAttenuation, linearAttenuation, quadraticAttenuation);
    }

    /// Returns the OpenGL index of this light source. range[0,MAX_LIGHTS]
    int getIndex(){
        return lightIndex_ - GL_LIGHT0;
    }

    //! Sets Curve object, which is used by followCurve().
    /// \param speed is being added to the current curve parameter each time
    ///           you call followCurve(). Must be inside interval ]0,1].
    /// \param curveParameter curve parameter followCurve() starts with.
    ///           Must be inside interval [0,1]
    /// \param start determines whether light's position is set to beginning of space curve
    void setCurve(Curve *curve, GLfloat speed = 0.01f, GLfloat curveParameter = 0.0f, bool start=true){
        curve_ = curve;

        if (start && curve != 0){
            setToCurvePoint(curveParameter);
        } else {
            setCurveParameter(curveParameter);
        }
        setSpeed(speed);

    }

    //! Get currently used Curve object.
    Curve* getCurve(){
        return curve_;
    }

    /// Sets the value, that is being added to the current curve parameter each time you call
    /// followCurve(). Must be inside interval ]0,1].
    void setSpeed(GLfloat speed){
        tgtAssert(speed > 0.f && speed <= 1.f, "Speed is expected to be inside the interval ]0,1]");
        speed_ = speed;
    }

    /// Gets current curve speed value.
    GLfloat getSpeed(){
        return speed_;
    }

    //! Sets current curve parameter.
    /// This value will be used the next time you call followCurve().
    /// Must be inside interval [0,1].
    void setCurveParameter(GLfloat curveParameter){
        tgtAssert(curveParameter >= 0.f && curveParameter <= 1.f,
            "Curve Parameter is expected to be inside the interval [0,1]");
        curveParameter_ = curveParameter;
    }

    /// Get current curve parameter.
    GLfloat getCurveParameter(){
        return curveParameter_;
    }


    //! Refreshes the light's position and spot direction.
    /// \attention {Call this method if the ModelView-Matrix
    ///                has changed and you want to reapply it to position and spot direction.}
    void refresh(){
        glLightfv(lightIndex_, GL_POSITION, position_.elem);
        glLightfv(lightIndex_, GL_SPOT_DIRECTION, spotDirection_.elem);
    }


    //! Sets light source's position to the point on the curve determined by current curve parameter
    /// and adds 'speed' to the curve parameter afterwards. \n
    /// If curve parameter exceeds 1.0, it is reset to 0.0.
    /// \note Of course, you have to set a curve first before you can use this method !
    void followCurve() {

        if (curve_ == 0) {
            tgtAssert(false, "No Curve specified.");
        } else {
            vec3 pos = curve_->getPoint(curveParameter_);
            vec4 homPos = vec4(pos, 1.f);
            setPosition(homPos);

            curveParameter_ += speed_;
            if (curveParameter_ > 1.f || curveParameter_ < 0.f){
                curveParameter_ = speed_;
            }
        }
    }

    //! Sets light source's position to the point on the curve determined by parameter t. \n
    /// The current curve parameter will also be set to t. \n
    /// Only values inside the range [0,1] are accepted.
    void setToCurvePoint(GLfloat t){
        setCurveParameter(t);
        followCurve();
    }


    /* Some administrative statics */

    /// Returns the maximum number of light sources supported by this OpenGL implementation.
    static int getMaxLights(){
        GLint maxLights;
        glGetIntegerv(GL_MAX_LIGHTS, &maxLights);
        return static_cast<int>(maxLights);
    }

    /// Returns the count of the currently instantiated Light objects.
    static int getLightCount(){
        return static_cast<int>(usedLightIndices_.size());
    }

    //! Returns true if OpenGL light limit is reached and no more Light objects can be created.
    /// Delete an existing Light object to be able to create a new one.
    static bool lightLimitReached(){
        return (getLightCount() >= getMaxLights());
    }



    /* Global light model switches */

    /// enable or disable lighting
    static void enableLighting(bool enable = true) {
        if (enable)
            glEnable(GL_LIGHTING);
        else
            glDisable(GL_LIGHTING);
    }

    /// is lighting enabled?
    static GLboolean isEnabledLighting() {
        return glIsEnabled(GL_LIGHTING);
    }

    /// sets a global ambient term for the scene
    static void setGlobalAmbient(const vec4& ambient) {
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient.elem);
    }

    /// gets the global ambient term of the scene
    static vec4 getGlobalAmbient() {
        vec4 ambient;
        glGetFloatv(GL_LIGHT_MODEL_AMBIENT, ambient.elem);
        return ambient;
    }

    //! enables or disables local viewer model \n
    /// See glLightModel documentation for an explanation.
    static void enableLocalViewer(bool enable = true) {
        glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, enable);
    }

    /// is local viewer model enabled?
    static GLboolean isEnabledLocalViewer() {
        GLboolean b;
        glGetBooleanv(GL_LIGHT_MODEL_LOCAL_VIEWER, &b);
        return b;
    }

    //! enables or disables two sided light model. \n
    /// if enabled, vertices of back face polygons are lighted using their back face materials and
    /// with reversed normals. \n
    /// See glLightModel documentation for additional information.
    static void enableTwoSide(bool enable = true) {
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, enable);
    }

    /// is two sided light model enabled?
    static GLboolean isEnabledTwoSide() {
        GLboolean b;
        glGetBooleanv(GL_LIGHT_MODEL_TWO_SIDE, &b);
        return b;
    }

    /**
     * sets the color control
     *
     * @param i either GL_SINGLE_COLOR or GL_SEPARATE_SPECULAR_COLOR
    **/
    static void setColorControl(GLint i) {
        tgtAssert( i == GL_SINGLE_COLOR || i == GL_SEPARATE_SPECULAR_COLOR,
            "Not a valid GL value for this operation" );
        glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, i);
    }

    /**
     * gets the color control
     *
     * @return either GL_SINGLE_COLOR or GL_SEPARATE_SPECULAR_COLOR
    **/
    static GLint getColorControl() {
        GLint i;
        glGetIntegerv(GL_LIGHT_MODEL_COLOR_CONTROL, &i);
        return i;
    }


    /// \cond Deprecated_Light

    /*
        DEPRECATED methods. Kept temporarily for downwards compatibility.
        Please change OpenGL light source parameters only by creating a Light object
        and calling its setter methods!!!
    */

    /// DEPRECATED method. Please change OpenGL light source parameters only via Light objects.
    static void enableLight(GLenum idx, bool enable = true) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        if (enable)
            glEnable(GL_LIGHT0 + idx);
        else
            glDisable(GL_LIGHT0 + idx);
    }

    /// DEPRECATED method. Please read OpenGL light source parameters only via Light objects.
    static GLboolean isEnabledLight(GLenum idx) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        return glIsEnabled(GL_LIGHT0 + idx);
    }

    /// DEPRECATED method. Please change OpenGL light source parameters only via Light objects.
    static void setPosition(GLenum idx, const vec4& position) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        glLightfv(GL_LIGHT0 + idx, GL_POSITION, position.elem);
    }

    /// DEPRECATED method. Please read OpenGL light source parameters only via Light objects.
    static vec4 getPosition(GLenum idx) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        vec4 position;
        glGetLightfv(GL_LIGHT0 + idx, GL_POSITION, position.elem);
        return position;
    }

    /// DEPRECATED method. Please change OpenGL light source parameters only via Light objects.
    static void setAmbient(GLenum idx, const vec4& ambient) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        glLightfv(GL_LIGHT0 + idx, GL_AMBIENT, ambient.elem);
    }

    /// DEPRECATED method. Please read OpenGL light source parameters only via Light objects.
    static vec4 getAmbient(GLenum idx) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        vec4 ambient;
        glGetLightfv(GL_LIGHT0 + idx, GL_AMBIENT, ambient.elem);
        return ambient;
    }

    /// DEPRECATED method. Please change OpenGL light source parameters only via Light objects.
    static void setDiffuse(GLenum idx, const vec4& diffuse) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        glLightfv(GL_LIGHT0 + idx, GL_DIFFUSE, diffuse.elem);
    }

    /// DEPRECATED method. Please read OpenGL light source parameters only via Light objects.
    static vec4 getDiffuse(GLenum idx) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        vec4 diffuse;
        glGetLightfv(GL_LIGHT0 + idx, GL_DIFFUSE, diffuse.elem);
        return diffuse;
    }

    /// DEPRECATED method. Please change OpenGL light source parameters only via Light objects.
    static void setSpecular(GLenum idx, const vec4& specular) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        glLightfv(GL_LIGHT0 + idx, GL_SPECULAR, specular.elem);
    }

    /// DEPRECATED method. Please read OpenGL light source parameters only via Light objects.
    static vec4 getSpecular(GLenum idx) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        vec4 specular;
        glGetLightfv(GL_LIGHT0 + idx, GL_SPECULAR, specular.elem);
        return specular;
    }

    /// DEPRECATED method. Please change OpenGL light source parameters only via Light objects.
    static void setConstantAttenuation(GLenum idx, float f) {
        glLightf(GL_LIGHT0 + idx, GL_CONSTANT_ATTENUATION, f);
    }

    /// DEPRECATED method. Please read OpenGL light source parameters only via Light objects.
    static float getConstantAttenuation(GLenum idx) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        float f;
        glGetLightfv(GL_LIGHT0 + idx, GL_CONSTANT_ATTENUATION, &f);
        return f;
    }

    /// DEPRECATED method. Please change OpenGL light source parameters only via Light objects.
    static void setLinearAttenuation(GLenum idx, float f) {
        glLightf(GL_LIGHT0 + idx, GL_LINEAR_ATTENUATION, f);
    }

    /// DEPRECATED method. Please read OpenGL light source parameters only via Light objects.
    static float getLinearAttenuation(GLenum idx) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        float f;
        glGetLightfv(GL_LIGHT0 + idx, GL_LINEAR_ATTENUATION, &f);
        return f;
    }

    /// DEPRECATED method. Please change OpenGL light source parameters only via Light objects.
    static void setQuadraticAttenuation(GLenum idx, float f) {
        glLightf(GL_LIGHT0 + idx, GL_QUADRATIC_ATTENUATION, f);
    }

    /// DEPRECATED method. Please read OpenGL light source parameters only via Light objects.
    static float getQuadraticAttenuation(GLenum idx) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        float f;
        glGetLightfv(GL_LIGHT0 + idx, GL_QUADRATIC_ATTENUATION, &f);
        return f;
    }

    /// DEPRECATED method. Please change OpenGL light source parameters only via Light objects.
    static void setSpotDirection(GLenum idx, const vec3& direction) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        glLightfv(GL_LIGHT0 + idx, GL_SPOT_DIRECTION, direction.elem);
    }

    /// DEPRECATED method. Please read OpenGL light source parameters only via Light objects.
    static vec3 getSpotDirection(GLenum idx) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        vec3 direction;
        glGetLightfv(GL_LIGHT0 + idx, GL_SPOT_DIRECTION, direction.elem);
        return direction;
    }

    /// DEPRECATED method. Please change OpenGL light source parameters only via Light objects.
    static void setSpotExponent(GLenum idx, float exponent) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        tgtAssert( exponent >= 0.f && exponent <= 128.f, "The exponent must be within [0, 128]" );
        glLightf(GL_LIGHT0 + idx, GL_SPOT_EXPONENT, exponent);
    }

    /// DEPRECATED method. Please read OpenGL light source parameters only via Light objects.
    static float getSpotExponent(GLenum idx) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        float exponent;
        glGetLightfv(GL_LIGHT0 + idx, GL_SPOT_EXPONENT, &exponent);
        return exponent;
    }

    /// DEPRECATED method. Please change OpenGL light source parameters only via Light objects.
    static void setSpotCutoff(GLenum idx, float cutoff) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        glLightf(GL_LIGHT0 + idx, GL_SPOT_CUTOFF, cutoff);
    }

   /// DEPRECATED method. Please read OpenGL light source parameters only via Light objects.
    static float getSpotCutoff(GLenum idx) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        float cutoff;
        glGetLightfv(GL_LIGHT0 + idx, GL_SPOT_CUTOFF, &cutoff);
        return cutoff;
    }


    /*
        create specific lights
    */

    /// DEPRECATED method. Please create Light object only by calling createLight().
    static void createDirectionalLight(GLenum idx, const vec3& direction) {
        tgtAssert( idx <= 7, "Wrong light source index" );

        // turn off spot lighting
        glLightf(GL_LIGHT0 + idx, GL_SPOT_CUTOFF, 180.f);

        // set a direction
        vec4 position(direction, 0.f);
        glLightfv(GL_LIGHT0 + idx, GL_POSITION, position.elem);
    }

    /// DEPRECATED method. Please create Light object only by calling createLight().
    static void createPointLight(GLenum idx, const vec3& position) {
        tgtAssert( idx <= 7, "Wrong light source index" );

        // turn off spot lighting
        glLightf(GL_LIGHT0 + idx, GL_SPOT_CUTOFF, 180.f);

        // set a postion
        vec4 pos(position, 1.f);
        glLightfv(GL_LIGHT0 + idx, GL_POSITION, position.elem);
    }

    /// DEPRECATED method. Please create Light object only by calling createLight().
    static void createSpotLight(GLenum idx, const vec3& position, const vec3& spotDirection, float cutoff, float exponent) {
        tgtAssert( idx <= 7, "Wrong light source index" );
        tgtAssert( (cutoff >= 0.f && exponent <= 90.f) || (cutoff == 180.f),
            "The cutoff value must be within [0, 90] or the special value 180" );
        tgtAssert( exponent >= 0.f && exponent <= 128.f, "The exponent must be within [0, 128]" );

        // set a postion
        vec4 pos(position, 1.f);
        glLightfv(GL_LIGHT0 + idx, GL_POSITION, pos.elem);

        // do the other spot light specific stuff
        glLightfv(GL_LIGHT0 + idx, GL_SPOT_DIRECTION, spotDirection.elem);
        glLightf (GL_LIGHT0 + idx, GL_SPOT_CUTOFF, cutoff);
        glLightf (GL_LIGHT0 + idx, GL_SPOT_EXPONENT, exponent);
    }

    /*
        End DEPRECATED methods.
    */

    /// \endcond


private:

    // index of the light source. used for glLight calls.
    int lightIndex_;

    // OpenGL light source attributes
    bool isEnabled_;
    vec4 position_;
    vec4 ambient_;
    vec4 diffuse_;
    vec4 specular_;
    vec3 spotDirection_;
    GLfloat spotExponent_;
    GLfloat spotCutoff_;
    GLfloat constantAttenuation_;
    GLfloat linearAttenuation_;
    GLfloat quadraticAttenuation_;

    // Curve object and parameters, used by followCurve()
    Curve *curve_;
    GLfloat curveParameter_;
    GLfloat speed_;

    // list of all light source indices currently in use.
    static std::list<int> usedLightIndices_;


    // private constructor.
    // called by createLight().
    Light(int index) {
        lightIndex_ = index;
    }

    // returns next currently unused light source index
    // or GL_LIGHT0-1, if there is no free index.
    // called by createLight().
    //
    //    Note: "It is always the case that GL_LIGHT$i$ = GL_LIGHT0 + $i$"
    //    says the OpenGL documentation.
    //
    static int getNextLightIndex(){

        int index = GL_LIGHT0-1;

        std::list<int>::iterator iter;
        for (int i = GL_LIGHT0; i < (GL_LIGHT0 + getMaxLights()) && index < GL_LIGHT0; i++ ){
            iter = find(usedLightIndices_.begin(), usedLightIndices_.end(), i);
            if (iter == usedLightIndices_.end()) index = i;
        }
        return index;
    }

};

} // namespace tgt

#endif //VRN_OPENGL_COMPATIBILITY_PROFILE

#endif //TGT_LIGHT_H
