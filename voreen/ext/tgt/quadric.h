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

#ifndef TGT_QUADRICS_H
#define TGT_QUADRICS_H

#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE

#include "tgt/tgt_gl.h"

#include "tgt/renderable.h"
#include "tgt/bounds.h"
#include "tgt/vector.h"
#include "tgt/light.h"
#include "tgt/material.h"

#include <math.h>

namespace tgt {

/*
    \todo Bounds!
*/

/**
    This class is a wrapper for GLU quadrics. It is
    meant to be a superclass for tgt quadrics and should not
    be instantiated. \n
    All tgt quadrics are inherited from this class and must
    implement the virtual render() method.
    \note { The overridden render() method should call Quadric's render()
            method at the beginning! }
    \note { Bounds are not implemented yet! }
*/
class Quadric : public Renderable {
public:

    /// This struct is used to specifie which kind of
    /// normals are generated.
    enum NormalStyle {
        NONE    = GLU_NONE,        /*!< No normals */
        FLAT    = GLU_FLAT,        /*!< One normal is generated for every facet of a quadric */
        SMOOTH  = GLU_SMOOTH    /*!< One normal is generated for every vertex of a quadric.
                                     This is the default value.*/
    };

    /// This struct is used to specifie the direction of the normals. ///
    enum Orientation {
        OUTSIDE = GLU_OUTSIDE ,    /*!< Draw quadric with normals pointing outward. This is the default value. */
        INSIDE  = GLU_INSIDE    /*!< Draw quadric with normals pointing inward. */
    };

    /// This struct is used to specifie the quadrics draw style.
    enum DrawStyle {
        FILL        = GLU_FILL,            /*!< Quadrics are rendered with polygon primitives.    */
        LINE        = GLU_LINE,            /*!< Quadrics are rendered as a set of lines.        */
        SILHOUETTE  = GLU_SILHOUETTE,    /*!< Quadrics are rendered as a set of lines,
                                             except that edges separating coplanar faces will not be drawn. */
        POINT       = GLU_POINT            /*!< Quadrics are rendered as a set of points.        */
    };

protected:

    GLUquadricObj*  quadric_;
    /// The quadric's material.
    Material* material_;
    bool materialEnabled_;

private:

    NormalStyle normalStyle_;        /*!< \see NormalStyle */
    bool        textureEnabled_;    /*!< specifies wether texture coordinates are generated */
    Orientation orientation_;        /*!< \see Orientation */
    DrawStyle   drawStyle_;            /*!< \see DrawStyle */

public:

    /// Constructor.
    Quadric(const Bounds& bounds, bool _static, bool visible)
      : Renderable      (bounds, _static, visible),
        quadric_        (gluNewQuadric()),
        material_       (NULL),
        materialEnabled_(false),
        normalStyle_    (SMOOTH),
        textureEnabled_ (false),
        orientation_    (OUTSIDE),
        drawStyle_      (FILL)
    {}

    virtual ~Quadric() {
        gluDeleteQuadric(quadric_);
    }

    virtual void render(){
        if (Light::isEnabledLighting() && materialEnabled_ && material_ != NULL){
            material_->activate();
        }
    }

    /* Getters and Setters */

    /// Return which kind of normals are generated.
    /// \see NormalStyle
    NormalStyle getNormalStyle() const {
        return normalStyle_;
    }

    /// Set which kind of normals are generated.
    /// \see NormalStyle
    void setNormalStyle(NormalStyle normalStyle) {
        gluQuadricNormals(quadric_, normalStyle);
        normalStyle_ = normalStyle;
    }

    /// Return whether quadric is to be textured.
    bool getTextureEnabled() const {
        return textureEnabled_;
    }

    /// Specifie whether quadric is to be textured.
    /// \param true Generate texture coordinates.
    /// \param false Do not generate texture coordinates. This is the default value.
    void setTextureEnabled(bool textureEnabled) {
        gluQuadricTexture(quadric_, textureEnabled);
        textureEnabled_ = textureEnabled;
    }

    /// Get orientation of normals.
    /// \see Orientation
    Orientation getOrientation() const {
        return orientation_;
    }

    /// Set orientation of normals.
    /// \see Orientation
    void setOrientation(Orientation orientation) {
        gluQuadricOrientation(quadric_, orientation);
        orientation_ = orientation;
    }

    /// Return the quadrics draw style.
    /// \see DrawStyle
    DrawStyle getDrawStyle() const {
        return drawStyle_;
    }

    /// Specifie the quadrics draw style.
    /// \see DrawStyle
    void setDrawStyle(DrawStyle drawStyle) {
        gluQuadricDrawStyle  (quadric_, drawStyle);
        drawStyle_ = drawStyle;
    }

    //! apply a material to the quadric
    /// \param enabled enable or disable use of material for rendering
    void setMaterial(Material* material, bool enabled = true){
        material_ = material;
        materialEnabled_ = enabled;
    }

    //! get currently used material
    Material* getMaterial(){
        return material_;
    }

    //! enable or disable use of material for rendering
    /// \note{ acitvation of material is only accepted, if material
    ///        has already been set }
    void enableMaterial(bool enable = true){
        if (enable){
            tgtAssert(material_ != NULL, "No material specified!");
            materialEnabled_ = true;
        } else {
            materialEnabled_ = false;
        }
    }

    //! check wether material is activated
    bool isEnabledMaterial(){
        return materialEnabled_;
    }

    virtual void init() {}
    virtual void deinit() {}
};

/**
    Is used to draw a sphere of the given radius centered around the origin.
*/
class Sphere : public Quadric {
protected:

    GLdouble    radius_;        /*!< radius of the sphere */
    GLint       slices_;        /*!< The number of subdivisions around the z-axis (similar to lines of longitude).    */
    GLint       stacks_;        /*!< The number of subdivisions along the z-axis (similar to lines of latitude).    */

public:

    /// Constructor.
    /// \param radius radius of the sphere
    /// \param slices The number of subdivisions around the z-axis (similar to lines of longitude).
    /// \param stacks The number of subdivisions along the z-axis (similar to lines of latitude).
    /// \param visible An invisible object is not rendered.
    /// \param _static A static object is not expected to move.
    Sphere( GLdouble radius, GLint slices, GLint stacks,
            bool visible = true, bool _static = true)

      : Quadric(Bounds(vec3(-(float)radius), vec3((float)radius)), _static, visible),
        radius_(radius),
        slices_(slices),
        stacks_(stacks)
    {}

    /// Renders the sphere around the origin.
    void render() {
        if (visible_) {
            Quadric::render();
            gluSphere(quadric_, radius_, slices_, stacks_);
        }
    }

    /// Set radius of the sphere.
    void setRadius(GLdouble radius){
        radius_ = radius;
    }

    /// Get radius of the sphere.
    GLdouble getRadius(){
        return radius_;
    }

    /// Set number of slices.
    /// \see slices_
    void setSlices(GLint slices){
        slices_ = slices;
    }

    /// Get number of slices.
    /// \see slices_
    GLint getSlices(){
        return slices_;
    }


    /// Set number of stacks.
    /// \see stacks_
    void setStacks(GLint stacks){
        stacks_ = stacks;
    }

    /// Get number of stacks.
    /// \see stacks_
    GLint getStacks(){
        return stacks_;
    }

};

/**
    A (two-dimensional) disk on the z = 0 plane
*/
class Disk : public Quadric {
protected:

    GLdouble    innerRadius_;    /*!< The inner radius of the disk (may be zero). If innerRadius > 0
                                     the disk has a hole of innerRadius in the middle.    */
    GLdouble    outerRadius_;    /*!< The outer radius of the disk.                        */
    GLint       slices_;        /*!< The number of subdivisions around the z-axis.        */
    GLint       loops_;            /*!< The number of concentric rings about the origin
                                     into which the disk is subdivided.                    */

public:

    /// Constructor.
    /// \param innerRadius \see innerRadius_
    /// \param outerRadius \see outerRadius_
    /// \param slices \see slices_
    /// \param loops \see loops_
    /// \param visible An invisible object is not rendered.
    /// \param _static A static object is not expected to move.
    Disk(   GLdouble innerRadius, GLdouble outerRadius, GLint slices, GLint loops,
            bool visible = true, bool _static = true)

      : Quadric(Bounds(), _static, visible),
        innerRadius_(innerRadius),
        outerRadius_(outerRadius),
        slices_(slices),
        loops_ (loops)
    {}

    /// Renders the disk on the z = 0 plane, centered at the origin.
    virtual void render() {
        if (visible_) {
            Quadric::render();
            gluDisk(quadric_, innerRadius_, outerRadius_, slices_, loops_);
        }
    }

    ///
    /// \see innerRadius_
    void setInnerRadius(GLdouble innerRadius){
        innerRadius_ = innerRadius;
    }

    /// \see innerRadius_
    GLdouble getInnerRadius(){
        return innerRadius_;
    }

    ///
    /// \see outerRadius_
    void setOuterRadius(GLdouble outerRadius){
        outerRadius_ = outerRadius;
    }

    ///
    /// \see outerRadius_
    GLdouble getOuterRadius(){
        return outerRadius_;
    }

    ///
    /// \see slices_
    void setSlices(GLint slices){
        slices_ = slices;
    }

    ///
    /// \see slices_
    GLint getSlices(){
        return slices_;
    }

    ///
    /// \see loops_
    void setLoops(GLint loops){
        loops_ = loops;
    }

    ///
    /// \see loops_
    GLint getLoops(){
        return loops_;
    }

};

/**
    A cylinder oriented along the z-axis. \n
    The base of the cylinder is placed at z = 0, and the top at z = height.
    By setting base or top radius to zero you get a cone.
    The cylinder is open: it has no bottom and top faces.
*/
class Cylinder : public Quadric {
protected:

    GLdouble    baseRadius_;    /*!< The radius of the cylinder at z = 0.            */
    GLdouble    topRadius_;        /*!< The radius of the cylinder at z = height_.        */
    GLdouble    height_;        /*!< The height of the cylinder.                    */
    GLint       slices_;        /*!< The number of subdivisions around the z-axis.    */
    GLint       stacks_;        /*!< The number of subdivisions along the z-axis.    */

public:

    /// Constructor.
    /// \param baseRadius \see baseRadius_
    /// \param topRadius \see topRadius_
    /// \param height \see height_
    /// \param slices \see slices_
    /// \param stacks \see stacks_
    /// \param visible An invisible object is not rendered.
    /// \param _static A static object is not expected to move.
      Cylinder(   GLdouble baseRadius, GLdouble topRadius, GLdouble height, GLint slices, GLint stacks,
                bool visible = true, bool _static = true)

      : Quadric(Bounds(), _static, visible),
        baseRadius_(baseRadius),
         topRadius_( topRadius),
        height_(height),
        slices_(slices),
        stacks_(stacks)
    {}

    /// Renders the cylinder along the z-axis.
    /// The base of the cylinder is placed at z = 0, and the top at z = height.
    virtual void render() {
        if (visible_) {
            Quadric::render();
            gluCylinder(quadric_, baseRadius_, topRadius_, height_, slices_, stacks_);
        }
    }

    ///
    /// \see baseRadius_
    void setBaseRadius(GLdouble baseRadius){
        baseRadius_ = baseRadius;
    }

    ///
    /// \see baseRadius_
    GLdouble getBaseRadius(){
        return baseRadius_;
    }

    ///
    /// \see topRadius_
    void setTopRadius(GLdouble topRadius){
        topRadius_ = topRadius;
    }

    ///
    /// \see topRadius_
    GLdouble getTopRadius(){
        return topRadius_;
    }

    ///
    /// \see height_
    void setHeight(GLdouble height){
        height_ = height;
    }

    ///
    /// \see height_
    GLdouble getHeight(){
        return height_;
    }

    ///
    /// \see slices_
    void setSlices(GLint slices){
        slices_ = slices;
    }

    ///
    /// \see slices_
    GLint getSlices(){
        return slices_;
    }

    ///
    /// \see stacks_
    void setStacks(GLint stacks){
        stacks_ = stacks;
    }

    ///
    /// \see stacks_
    GLint getStacks(){
        return stacks_;
    }
};

/**
    This is the solid version of Cylinder: it is rendered with bottom and top faces.
    The base of the cylinder is placed at z = 0, and the top at z = height.
    By setting base or top radius to zero you get a solid cone.
*/
class SolidCylinder : public Cylinder {
protected:

    Disk bottom_;
    Disk top_;

public:

    /// Constructor.
    /// \param baseRadius The radius of the cylinder at z = 0.
    /// \param topRadius The radius of the cylinder at z = height.
    /// \param height Height of the cylinder.
    /// \param slices The number of subdivisions around the z-axis. Applies to the cylinder and the
    ///                    bottom and top faces.
    /// \param stacks The number of subdivisions of the cylinder along the z-axis.
    /// \param loops The number of concentric rings about the origin into which the bottom and top disks are subdivided.
    /// \param visible An invisible object is not rendered.
    /// \param _static A static object is not expected to move.
    SolidCylinder(   GLdouble baseRadius, GLdouble topRadius, GLdouble height, GLint slices, GLint stacks,
                GLint loops, bool visible = true, bool _static = true)

        : Cylinder(baseRadius, topRadius, height, slices, stacks, _static, visible),
          bottom_(0, baseRadius, slices, loops, _static, visible),
          top_(0, topRadius, slices, loops, _static, visible)
    {
        bottom_.setOrientation(Quadric::INSIDE);
    }

    /// Renders the cylinder along the z-axis.
    /// The base of the cylinder is placed at z = 0, and the top at z = height.
    virtual void render(){
        if (visible_){
            Quadric::render();
            Cylinder::render();
            if (baseRadius_ > 0){
                bottom_.render();
            }
            if (topRadius_ > 0){
                glMatrixMode(GL_MODELVIEW);
                glPushMatrix();
                glTranslated(0,0,height_);
                top_.render();
                glPopMatrix();
            }
        }
    }

    void setBaseRadius(GLdouble baseRadius){
        baseRadius_ = baseRadius;
        bottom_.setOuterRadius(baseRadius);
    }

    void setTopRadius(GLdouble topRadius){
        topRadius_ = topRadius;
        top_.setOuterRadius(topRadius);
    }

    void setSlices(GLint slices){
        slices_ = slices;
        bottom_.setSlices(slices);
        top_.setSlices(slices);
    }

    void setLoops(GLint loops){
        bottom_.setLoops(loops);
        top_.setLoops(loops);
    }

    GLint getLoops(){
        return bottom_.getLoops();
    }

    void setStatic(bool _static){
        Cylinder::setStatic(_static);
        bottom_.setStatic(_static);
        top_.setStatic(_static);
    }

    void setNormalStyle(NormalStyle normalStyle) {
        Cylinder::setNormalStyle(normalStyle);
        bottom_.setNormalStyle(normalStyle);
        top_.setNormalStyle(normalStyle);
    }

    void setTextureEnabled(bool _enabled){
        Cylinder::setTextureEnabled(_enabled);
        bottom_.setTextureEnabled(_enabled);
        top_.setTextureEnabled(_enabled);
    }

    void setOrientation(Orientation orientation){
        Cylinder::setOrientation(orientation);
        if (orientation == OUTSIDE)
            bottom_.setOrientation(INSIDE);
        else
            bottom_.setOrientation(OUTSIDE);
        top_.setOrientation(orientation);
    }

    void setDrawStyle(DrawStyle drawStyle){
        Cylinder::setDrawStyle(drawStyle);
        bottom_.setDrawStyle(drawStyle);
        top_.setDrawStyle(drawStyle);
    }
};


/**
    A rectangle oriented at x-y-plane. \n
    The center of the rectangle is placed at x = 0, y = 0, z = 0.
*/
class Rect : public Quadric {
protected:

    GLdouble    width_;
    GLdouble    height_;

    GLdouble    h_width_;
    GLdouble    h_height_;

public:

    /// Constructor.
    /// \param width The width of the rectangle.
    /// \param height The height of the rectangle.
    /// \param visible An invisible object is not rendered.
    /// \param _static A static object is not expected to move.
    Rect(   GLdouble width, GLdouble height,
        bool visible = true, bool _static = true)
      : Quadric(Bounds(), _static, visible)
    {
        setWidth(width);
        setHeight(height);
    }

    /// Renders the rectangle in the x-y-plane.
    /// The center of the rectangle is in (0,0,0)
    /// e.g. the measures in x direction: -(width / 2.0), (width / 2.0)
    virtual void render(){
        if (visible_) {

            Quadric::render();

            switch (getDrawStyle()) {
            case LINE:
                glBegin(GL_LINE_STRIP);
                break;
            case SILHOUETTE:
                glBegin(GL_LINE_STRIP);
                break;
            case POINT:
                glBegin(GL_POINTS);
                break;
            default:
                glBegin(GL_QUADS);
                break;
            }


            // if FLAT or SMOOTH set normals
            if (getNormalStyle() != NONE) {
                //set normal (inside or outside)
                if (getOrientation() == INSIDE) {
                    glNormal3f(0,0,-1);
                } else {
                    glNormal3f(0,0, 1);
                }
            }

            // render the rectangle
            if (getTextureEnabled())
                glTexCoord2f(0.f, 0.f);
            glVertex3d(-h_width_, -h_height_, 0);
            if (getTextureEnabled())
                glTexCoord2f(0.f, 1.f);
            glVertex3d(-h_width_,  h_height_, 0);
            if (getTextureEnabled())
                glTexCoord2f(1.f, 1.f);
            glVertex3d( h_width_,  h_height_, 0);
            if (getTextureEnabled())
                glTexCoord2f(1.f, 0.f);
            glVertex3d( h_width_, -h_height_, 0);

            if ((getDrawStyle() == LINE) || (getDrawStyle() == SILHOUETTE)) {
                glVertex3d(-h_width_, -h_height_, 0);
            }

            glEnd();
        }
    }

    /// Returns the width of the rectangle
    GLdouble getWidth() {
        return width_;
    }

    /// Sets the width of the rectangle
    void setWidth(GLdouble w) {
        width_ = w;
        h_width_ = w / 2.0;
    }

    /// Returns the height of the rectangle
    GLdouble getHeight() {
        return height_;
    }

    /// Sets the height of the rectangle
    void setHeight(GLdouble h) {
        height_ = h;
        h_height_ = h / 2.0;
    }
};

/**
    A quad centered around the origin. \n
    It is created of six rectangle objects.
*/
class Quad : public Quadric {
protected:

    GLdouble    width_;
    GLdouble    height_;
    GLdouble    depth_;

private:

    GLdouble    h_width_;
    GLdouble    h_height_;
    GLdouble    h_depth_;

    Rect* front_;
    Rect* back_;
    Rect* top_;
    Rect* bottom_;
    Rect* right_;
    Rect* left_;

    void initializeRects() {
            front_  = new Rect(width_, height_, true);
            back_   = new Rect(width_, height_, true);
            top_    = new Rect(width_, depth_, true);
            bottom_ = new Rect(width_, depth_, true);
            right_  = new Rect(height_, depth_, true);
            left_   = new Rect(height_, depth_, true);
    }

    void deleteRects() {
            delete front_;
            delete back_;
            delete top_;
            delete bottom_;
            delete right_;
            delete left_;
    }

public:

    /// Constructor.
    /// \param width The width of the rectangle.
    /// \param height The height of the rectangle.
    /// \param depth The depth of the rectangle.
    /// \param visible An invisible object is not rendered.
    /// \param _static A static object is not expected to move.
    Quad(   GLdouble width, GLdouble height, GLdouble depth,
        bool visible = true, bool _static = true)
      : Quadric(Bounds(), _static, visible)
    {
            initializeRects();

            setWidth(width);
            setHeight(height);
            setDepth(depth);
    }

        /// Destructor
        ~Quad() {
            deleteRects();
        }

    /// Renders the quad.
    /// The center of the quad is in (0,0,0)
    /// eg. the measures in x direction: -(width / 2.0), (width / 2.0)
    virtual void render(){
        if (visible_) {

            Quadric::render();

            glPushMatrix();
                glTranslated(0,0,h_depth_);
                front_->render();
            glPopMatrix();

            glPushMatrix();
                glTranslated(0,0,-h_depth_);
                glRotatef(180.f, 1.f, 0.f, 0.f);
                back_->render();
            glPopMatrix();

            glPushMatrix();
                glTranslated(0,h_height_,0);
                glRotatef(90.f, -1.f, 0.f, 0.f);
                top_->render();
            glPopMatrix();

            glPushMatrix();
                glTranslated(0,-h_height_,0);
                glRotatef(90.f, 1.f, 0.f, 0.f);
                bottom_->render();
            glPopMatrix();

            glPushMatrix();
                glTranslated(h_width_,0,0);
                glRotatef(90.f, 0.f, 1.f, 0.f);
                glRotatef(90.f, 0.f, 0.f, 1.f);
                right_->render();
            glPopMatrix();

            glPushMatrix();
                glTranslated(-h_width_,0,0);
                glRotatef(90, 0.f,-1.f, 0.f);
                glRotatef(90, 0.f,0.f, 1.f);
                left_->render();
            glPopMatrix();
        }
    }

    /// Returns the width of the quad
    GLdouble getWidth() {
        return width_;
    }

    /// Sets the width of the quad
    void setWidth(GLdouble w) {
        width_ = w;
        h_width_ = width_ / 2.0;
                deleteRects();
        initializeRects();
    }

    /// Returns the height of the quad
    GLdouble getHeight() {
        return height_;
    }

    /// Sets the height of the quad
    void setHeight(GLdouble h) {
        height_ = h;
        h_height_ = height_ / 2.0;
                deleteRects();
        initializeRects();
    }

    /// Returns the depth of the quad
    GLdouble getDepth() {
        return depth_;
    }

    /// Sets the depth of the quad
    void setDepth(GLdouble d) {
        depth_ = d;
        h_depth_ = depth_ / 2.0;
                deleteRects();
        initializeRects();
    }

    NormalStyle getNormalStyle() {
        return front_->getNormalStyle();
    }

    void setNormalStyle(NormalStyle normalStyle) {
        front_->setNormalStyle(normalStyle);
        back_->setNormalStyle(normalStyle);
        bottom_->setNormalStyle(normalStyle);
        top_->setNormalStyle(normalStyle);
        right_->setNormalStyle(normalStyle);
        left_->setNormalStyle(normalStyle);
    }

    bool getTextureEnabled() {
        return front_->getTextureEnabled();
    }

    void setTextureEnabled(bool _enabled){
        front_->setTextureEnabled(_enabled);
        back_->setTextureEnabled(_enabled);
        bottom_->setTextureEnabled(_enabled);
        top_->setTextureEnabled(_enabled);
        right_->setTextureEnabled(_enabled);
        left_->setTextureEnabled(_enabled);
    }

    Orientation getOrientation() {
        return front_->getOrientation();
    }

    void setOrientation(Orientation orientation){
        front_->setOrientation(orientation);
        back_->setOrientation(orientation);
        bottom_->setOrientation(orientation);
        top_->setOrientation(orientation);
        right_->setOrientation(orientation);
        left_->setOrientation(orientation);
    }

    DrawStyle getDrawStyle() {
        return front_->getDrawStyle();
    }

    void setDrawStyle(DrawStyle drawStyle){
        front_->setDrawStyle(drawStyle);
        back_->setDrawStyle(drawStyle);
        bottom_->setDrawStyle(drawStyle);
        top_->setDrawStyle(drawStyle);
        right_->setDrawStyle(drawStyle);
        left_->setDrawStyle(drawStyle);
    }
};


/**
    An arc of a disk specified by a start and a sweep angle. \n
    It is rendered on the z=0 plane, centered at the origin.
*/
class PartialDisk : public Disk {
protected:

    GLdouble    startAngle_;    /*!< The starting angle, in degrees, of the disk portion.
                                     0 degrees is along the positive y-axis.                */
    GLdouble    sweepAngle_;    /*!< The sweep angle, in degrees, of the disk portion.
                                     Rotation is clockwise around the z-axis.                */

public:

    /// Constructor.
    /// \param startAngle The starting angle, in degrees, of the disk portion. 0 degrees is along the positive y-axis.
    /// \param sweepAngle The sweep angle, in degrees, of the disk portion. Rotation is clockwise around the z-axis.
    PartialDisk(    GLdouble innerRadius, GLdouble outerRadius, GLint slices, GLint loops,
                    GLdouble startAngle, GLdouble sweepAngle,
                    bool visible = true, bool _static = true)
      : Disk(   innerRadius, outerRadius, slices, loops,
                _static, visible),
        startAngle_(startAngle),
        sweepAngle_(sweepAngle)
    {}

    /// Renders the partial disk on the z=0 plane, centered at the origin.
    virtual void render() {
        if (visible_){
            Quadric::render();
            gluPartialDisk( quadric_, innerRadius_, outerRadius_, slices_, loops_,
                            startAngle_, sweepAngle_);
        }
    }

    ///
    /// see \startAngle_
    void setStartAngle(GLdouble startAngle){
        startAngle_ = startAngle;
    }

    ///
    /// see \startAngle_
    GLdouble getStartAngle(){
        return startAngle_;
    }

    ///
    /// see \sweepAngle_
    void setSweepAngle(GLdouble sweepAngle){
        sweepAngle_ = sweepAngle;
    }

    ///
    /// see \sweepAngle_
    GLdouble getSweepAngle(){
        return sweepAngle_;
    }

};

/**
    An equilateral traingle oriented at x-y-plane. \n
    The center (i.e. center of inscribed circle) of the rectangle is placed at x = 0, y = 0, z = 0.
*/
class Triangle : public Quadric {
protected:

    GLdouble    size_;

public:

    /// Constructor.
    /// \param size The size of the triangle.
    /// \param visible An invisible object is not rendered.
    /// \param _static A static object is not expected to move.
    Triangle(   GLdouble size, bool visible = true, bool _static = true)
      : Quadric(Bounds(), _static, visible)
    {
        setSize(size);
    }

    /// Renders the triangle in the x-y-plane.
    /// The center of the rectangle is in (0,0,0)
    virtual void render(){
        if (visible_) {

            Quadric::render();

            glPushMatrix();
            //move the center to the right position
            glTranslated(0,-size_*sqrt(3.0)/6,0);
            switch (getDrawStyle()) {
            case LINE:
                glBegin(GL_LINE_STRIP);
                break;
            case SILHOUETTE:
                glBegin(GL_LINE_STRIP);
                break;
            case POINT:
                glBegin(GL_POINTS);
                break;
            default:
                glBegin(GL_TRIANGLES);
                break;
            }


            // if FLAT or SMOOTH set normals
            if (getNormalStyle() != NONE) {
                //set normal (inside or outside)
                if (getOrientation() == INSIDE) {
                    glNormal3f(0,0,-1);
                } else {
                    glNormal3f(0,0, 1);
                }
            }

            // render the triangle
            if (getTextureEnabled())
                glTexCoord2f(1.f, 0.f);
            glVertex3d(size_/2,0,0);
            if (getTextureEnabled())
                glTexCoord2f(0.5f, 1.f);
            glVertex3d(0,size_*sqrt(3.0)/2,0);
            if (getTextureEnabled())
                glTexCoord2f(0.f, 0.f);
            glVertex3d(-size_/2,0,0);

            if ((getDrawStyle() == LINE) || (getDrawStyle() == SILHOUETTE)) {
                glVertex3d(size_/2,0,0);
            }

            glEnd();
            glPopMatrix();
        }
    }

    /// Returns the size of the triangle
    GLdouble getSize() {
        return size_;
    }

    /// Sets the size of the triangle
    void setSize(GLdouble s) {
        size_ = s;
    }

};

/**
    A tetrahedron centered around the origin. \n
    It is created of four triangle objects.
*/
class Tetrahedron : public Quadric {
protected:

    GLdouble size_;

private:

    Triangle* front_;
    Triangle* bottom_;
    Triangle* right_;
    Triangle* left_;

    void initializeTriangles() {
        front_  = new Triangle(size_);
        bottom_ = new Triangle(size_);
        right_  = new Triangle(size_);
        left_   = new Triangle(size_);
    }

    void deleteTriangles() {
        delete front_;
        delete bottom_;
        delete right_;
        delete left_;
    }

public:

    /// Constructor.
    /// \param size The size of the tetrahedron.
    /// \param visible An invisible object is not rendered.
    /// \param _static A static object is not expected to move.
    Tetrahedron(   GLdouble size, bool visible = true, bool _static = true)
      : Quadric(Bounds(), _static, visible)
    {
            initializeTriangles();
            setSize(size);
    }

        /// Destructor
        ~Tetrahedron() {
            deleteTriangles();
        }

    /// Renders the etrahedron.
    /// The center of the etrahedron is in (0,0,0)
    virtual void render(){
        if (visible_) {

            Quadric::render();
            // for angles, lengths etc. see http://en.wikipedia.org/wiki/Tetrahedron
            float faceAngle = 70.5288;
            double length = size_*sqrt(3.0)/6;

            glPushMatrix();
            glTranslated(0,0,-size_/sqrt(24.0));

            glPushMatrix();
                glTranslated(0,-length,0);
                glRotatef(faceAngle,1.0,0,0);
                glTranslated(0,length,0);
                front_->render();
            glPopMatrix();

            glPushMatrix();
                bottom_->render();
            glPopMatrix();

            glPushMatrix();
                glRotatef(120,0,0,1);
                glTranslated(0,-length,0);
                glRotatef(faceAngle,1.0,0,0);
                glTranslated(0,length,0);
                right_->render();
            glPopMatrix();

            glPushMatrix();
                glRotatef(240,0,0,1);
                glTranslated(0,-length,0);
                glRotatef(faceAngle,1.0,0,0);
                glTranslated(0,length,0);
                left_->render();
            glPopMatrix();

            glPopMatrix();
        }
    }

    /// Returns the size of the tetrahedron
    GLdouble getSize() {
        return size_;
    }

    /// Sets the size of the tetrahedron
    void setSize(GLdouble s) {
        size_ = s;
        deleteTriangles();
        initializeTriangles();
    }

    NormalStyle getNormalStyle() {
        return front_->getNormalStyle();
    }

    void setNormalStyle(NormalStyle normalStyle) {
        front_->setNormalStyle(normalStyle);
        bottom_->setNormalStyle(normalStyle);
        right_->setNormalStyle(normalStyle);
        left_->setNormalStyle(normalStyle);
    }

    bool getTextureEnabled() {
        return front_->getTextureEnabled();
    }

    void setTextureEnabled(bool _enabled){
        front_->setTextureEnabled(_enabled);
        bottom_->setTextureEnabled(_enabled);
        right_->setTextureEnabled(_enabled);
        left_->setTextureEnabled(_enabled);
    }

    Orientation getOrientation() {
        return front_->getOrientation();
    }

    void setOrientation(Orientation orientation){
        front_->setOrientation(orientation);
        bottom_->setOrientation(orientation);
        right_->setOrientation(orientation);
        left_->setOrientation(orientation);
    }

    DrawStyle getDrawStyle() {
        return front_->getDrawStyle();
    }

    void setDrawStyle(DrawStyle drawStyle){
        front_->setDrawStyle(drawStyle);
        bottom_->setDrawStyle(drawStyle);
        right_->setDrawStyle(drawStyle);
        left_->setDrawStyle(drawStyle);
    }
};


}

#endif //VRN_OPENGL_COMPATIBILITY_PROFILE

#endif //TGT_QUADRICS_H
