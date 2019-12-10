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

#ifndef VRN_LABELINGMATH_H
#define VRN_LABELINGMATH_H

#include <cmath>
#include <vector>
#include "tgt/tgt_math.h"
#include "tgt/vector.h"
#include "tgt/matrix.h"
#include "tgt/tgt_gl.h"
#include "tgt/tgt_math.h"

namespace labeling {

/// Specifies a polynomial by its degree and
/// coefficients:
/// poly(x) = coeff[0] + coeff[1]*x + ...
struct Polynomial {
    int degree;
    float* coeff;
};

//----------------------------------------------------------------------------

/**
 * Superclass for 3-dimensional space curves.
 * f(t) = ( x(t), y(t), z(t) )
 * The curve is determined by an array of control points.
 */
class Curve3D {
public:

    /**
     * Default constructor.
     */
    Curve3D();

    /**
     * Frees all allocated resources.
     */
    virtual ~Curve3D() {}

    /**
     * Sets the control points of the curve. Internally the
     * curve's parameters are calculated from these control points.
     *
     * @param curveLength the length of the curve segment, whose parameter range
     *      range is mapped to [0;1]. The segment is centered around the control points.
     *      If no length is supplied, the segment runs from the first to the last
     *      control point (these points are only approximated, since this is a
     *      least-squares-fit)
     * @return true, if parameter calculation was successful
     */
    bool setCtrlPoints(std::vector<tgt::vec3> &ctrlPoints, float curveLength=0.f);

    /// Returns the control point array
    std::vector<tgt::vec3> getCtrlPoints(){
        return ctrlPoints_;
    }

    /// Returns the curve point at parameter \c t.
    virtual tgt::vec3 getCurvePoint(float t)=0;

    /// Returns the tangent vector at parameter \c t.
    virtual tgt::vec3 getTangent(float t)=0;

    /// Returns the length of the tangent vector at parameter \c t.
    virtual float getTangentMagnitude(float t)=0;

    /**
     *  Returns the length of the curve segment between
     *  parameters t1 and t2.
     */
    virtual float getSegmentLength(float t1, float t2)=0;

    /** Returns the curve point at parameter t that satisfies the following condition:
     *  - getSegmentLength(curParam, t) == offset
     * \note This is just an approximation.
     *
     * @param c curParam is set to t afterwards.
     */
    virtual tgt::vec3 getNextPoint(float &curParam, float offset)=0;

    /**
     * Shifts the whole curve by \c shiftVector:
     */
    virtual void shift(tgt::vec3 shiftVector)=0;

protected:

    /// Calculates the curve from the prior set control points.
    /// @param curveLength specifies the length of the curve segment whose
    ///     parameter range is mapped to [0;!]
    virtual bool calcFunction(float curveLength=0)=0;

    /// The Curve's control points.
    std::vector<tgt::vec3> ctrlPoints_;

    /// Scale and shift parameters for mapping parameter range of curve segment
    /// to [0;1]
    float paramScale_, paramShift_;

};

//----------------------------------------------------------------------------

/**
 * Superclass for 2-dimensional space curves: f(t) = ( x(t), y(t) ).
 * The curve is determined by an array of control points.
 */
class Curve2D {
public:

    /// Default constructor.
    Curve2D();

    /// Destructor. Frees all allocated resources.
    virtual ~Curve2D() {}

    /**
     * Sets the control points of the curve. Internally the
     * curve's parameters are calculated from these control points.
     *
     * @param curveLength the length of the curve segment, whose parameter range
     *      range is mapped to [0;1]. The segment is centered around the control points.
     *      If no length is supplied, the segment runs from the first to the last
     *      control point (these points are only approximated, since this is a
     *      least-squares-fit)
     * @return true, if parameter calculation was successful
     */
    bool setCtrlPoints(std::vector<tgt::vec2> &ctrlPoints, float curveLength=0.f);

    /// Returns the control point array
    std::vector<tgt::vec2> getCtrlPoints(){
        return ctrlPoints_;
    }

    /// Returns the curve point at parameter t
    virtual tgt::vec2 getCurvePoint(float t)=0;

    /// Returns the tangent vector at parameter t
    virtual tgt::vec2 getTangent(float t)=0;

    /// Returns the magnitude of the tangent vector at parameter t
    virtual float getTangentMagnitude(float t)=0;

    /**
     *  Returns the length of the curve segment between
     *  parameters t1 and t2.
     */
    virtual float getSegmentLength(float t1, float t2)=0;

    /** Returns the curve point at parameter t that satisfies the following condition:
     *  - getSegmentLength(curParam, t) == offset
     * \note This is just an approximation.
     *
     * @param c curParam is set to t afterwards.
     */
    virtual tgt::vec2 getNextPoint(float &curParam, float offset)=0;

    /// Shifts the whole curve by shiftVector.
    virtual void shift(tgt::vec2 shiftVector)=0;

protected:

    /// Calculates the curve from the prior set control points.
    /// @param curveLength specifies the length of the curve segment whose
    ///     parameter range is mapped to [0;!]
    virtual bool calcFunction(float curveLength=0)=0;

    /// The curve's control points
    std::vector<tgt::vec2> ctrlPoints_;

    /// Scale and shift parameters for mapping parameter range of curve segment
    /// to [0;1]
    float paramScale_, paramShift_;

};

//----------------------------------------------------------------------------

/**
 * Implementation of Curve3D.
 * Curve is a polynomial least-squares-fit to the control points.
 */
class Curve3DPolynomial : public Curve3D {
public:

    /// Creates the 3D polynomial.
    /// @param degree degree of polynomial
    Curve3DPolynomial(int degree);

    /// Destructor. Frees all allocated resources.
    virtual ~Curve3DPolynomial() {}

    tgt::vec3 getCurvePoint(float t);
    tgt::vec3 getTangent(float t);
    float getTangentMagnitude(float t);
    float getSegmentLength(float t1, float t2);
    tgt::vec3 getNextPoint(float &curParam, float offset);

    /// Shifts the whole curve by shiftVector.
    void shift(tgt::vec3 shiftVector);

    /// Returns coefficients of x(t) polynomial.
    float* getCoefficientsX(){
        return xfunc_.coeff;
    };

    /// Returns coefficients of y(t) polynomial.
    float* getCoefficientsY(){
        return yfunc_.coeff;
    };

    /// Returns coefficients of z(t) polynomial.
    float* getCoefficientsZ(){
        return zfunc_.coeff;
    };

protected:

    bool calcFunction(float curveLength=0);

    int degree_;
    int numCoeff_;
    Polynomial xfunc_;
    Polynomial yfunc_;
    Polynomial zfunc_;

};

//----------------------------------------------------------------------------

/**
 * Implementation of Curve2D.
 * Curve is a polynomial least-squares-fit to the control points.
 */
class Curve2DPolynomial : public Curve2D {
public:

    /// Creates the 3D polynomial
    /// @param degree degree of polynomial
    Curve2DPolynomial(int degree);

    /// Destructor. Frees all allocated resources.
    virtual ~Curve2DPolynomial() {}

    tgt::vec2 getCurvePoint(float t);
    tgt::vec2 getTangent(float t);
    float getTangentMagnitude(float t);
    float getSegmentLength(float t1, float t2);
    tgt::vec2 getNextPoint(float &curParam, float offset);

    /// Shifts the whole curve by shiftVector.
    void shift(tgt::vec2 shiftVector);

    /// Returns coefficients of x(t) polynomial.
    float* getCoefficientsX(){
        return xfunc_.coeff;
    };

    /// Returns coefficients of y(t) polynomial.
    float* getCoefficientsY(){
        return yfunc_.coeff;
    };

private:

    bool calcFunction(float curveLength=0);

    int degree_;
    int numCoeff_;
    Polynomial xfunc_;
    Polynomial yfunc_;

};

//----------------------------------------------------------------------------

/// A bezier surface with arbitrary dimensions.
/// The patch is determined by an array of control points with dimensions
/// (degreeS+1)x(degreeT+1).
/// The surface point at parameters s=t=0 is the lower left control point,
/// the surface point at parameters s=t=1 is the upper right control point.
class BezierPatch {
public:

    /// Default constructor.
    BezierPatch(bool useDisplayList=true);

    /// Default destructor. Frees all allocated resources.
    ~BezierPatch();

    /// Sets the patch's control points.
    /// @param ctrlPoints control point array, ordered row-wise, beginning at the
    ///     bottom row
    /// @param degreeS patch's degree in "horizontal" direction:
    ///     length of a row of the control point array = degreeS+1
    /// @param degreeT patch's degree in "vertical" direction:
    ///     length of a column of the control point array = degreeT+1
    void setCtrlPoints(tgt::vec3* ctrlPoints, int degreeS, int degreeT);

    /// Returns the patch's horizontal dimension (see constructor).
    int getDegreeS();

    /// Returns the patch's vertical dimension (see constructor)
    int getDegreeT();

    /// Returns a pointer to the control point array and writes the patch's
    /// dimensions to ref parameters degreeS and degreeT.
    tgt::vec3* getCtrlPoints(int &degreeS, int &degreeT);

    /// Returns the surface point at parameters s and t.
    tgt::vec3 getPoint(float s, float t);

    /// Returns the tangent in s-direction at the surface point
    /// at parameters s and t.
    tgt::vec3 getTangentS(float s, float t);

    /// Returns the tangent in t-direction at the surface point
    /// at parameters s and t.
    tgt::vec3 getTangentT(float s, float t);

    /// Returns the normal at the surface point
    /// at parameters s and t.
    tgt::vec3 getNormal(float s, float t);

    /// Renders the whole patch (parameter range [0:1]^2) as a quad strip
    /// and optionally generates texture coordinates for an arbitrary texture unit.
    /// Internally a display list is used to speed up rendering.
    /// @param s_steps number of samplings in s-direction (horizontal)
    /// @param t_steps number of samplings in t-direction (vertical)
    /// @param genTextureCoords determines if texture coordinates are generated
    /// @param texUnit texture unit the texture coordinates are generated for
    void render(int s_steps = 10, int t_steps = 10,
        bool genTexCoords=false, GLuint texUnit=0);

private:
    int degreeS_;
    int degreeT_;
    tgt::vec3* ctrlPoints_;

    GLuint displayList_;
    bool useDisplayList_;
};


//----------------------------------------------------------------------------

/// Template class representing a bitmap.
/// Saves bitmap data and provides functions for accessing and manipulating
/// it on a per-pixel/per-subpixel basis.
/// Note:  For performance reasons none of the functions do parameter checking.
///        Obeying the specifications is up to the function's callers!
template<class T>
class Bitmap {
public:

    /// Creates an empty bitmap with dimension 0x0.
    Bitmap() {
        data_ = NULL;
        width_ = 0;
        height_ = 0;
        channels_ = 1;
    };

    /// Creates an empty bitmap with dimensions width x height and
    /// the specified number of channels. The necessary memory
    /// is allocated.
    Bitmap(int width, int height, int channels = 1) {
        data_ = new T[width*height*channels];
        width_ = width;
        height_ = height;
        channels_ = channels;
    };

    /// Creates a bitmap by providing the data as an array and specifying
    /// its dimensions and channels count.
    Bitmap(T* data, int width, int height, int channels = 1) {
        data_ = data;
        width_ = width;
        height_ = height;
        channels_ = channels;
    };

    ~Bitmap(){
        delete[] data_;
        data_ = 0;
    }

    /// Sets the bitmap's data by providing the data array and
    /// its dimensions and channels count.
    void setData(T* data, int width, int height, int channels = 1){
        if(data_) {
            delete[] data_;
            data_ = 0;
        }
        data_ = data;
        width_ = width;
        height_ = height;
        channels_ = channels;
    }

    /// Returns a pointer to the data.
    T* getData(){
        return data_;
    }

    /// Sets the element at pixel coordinates x and y and
    /// the specified channel.
    inline void setElem(int x, int y, T elem, int channel=0){
        setElem((width_*y + x)*channels_ + channel, elem);
    }

    /// Sets the element at array position index.
    /// Use not recommended!
    inline void setElem(int index, T elem){
        data_[index] = elem;
    }

    /// Returns the element at pixel coordinates x and y and
    /// the specified channel.
    inline T getElem(int x, int y, int channel=0){
        return getElem((width_*y + x)*channels_ + channel);
    }

    /// Returns the element at array position index.
    /// Use not recommended!
    inline T getElem(int index){
        return data_[index];
    }

    /// Returns a 2D-Vector consisting of the first 2 channels
    /// of pixel (x,y).
    inline tgt::Vector2<T> getPixel2Ch(int x, int y){
        return tgt::Vector2<T>(data_[(width_*y + x)*channels_],
            data_[(width_*y + x)*channels_ + 1] );
    }

    /// Returns a 3D-Vector consisting of the first 3 channels
    /// of pixel (x,y).
    inline tgt::Vector3<T> getPixel3Ch(int x, int y){
        x = tgt::clamp(x, 0, width_-1);
        y = tgt::clamp(y, 0, height_-1);
        return tgt::Vector3<T>(data_[(width_*y + x)*channels_],
            data_[(width_*y + x)*channels_ + 1],
            data_[(width_*y + x)*channels_ + 2]  );
    }

    /// Returns a 4D-Vector consisting of the first 4 channels
    /// of pixel (x,y).
    inline tgt::Vector4<T> getPixel4Ch(int x, int y){
        return tgt::Vector4<T>(data_[(width_*y + x)*channels_],
            data_[(width_*y + x)*channels_ + 1],
            data_[(width_*y + x)*channels_ + 2],
            data_[(width_*y + x)*channels_ + 3] );
    }

    /// Merges this bitmap with src: This bitmap's pixels are set to src's pixels.
    /// If src's dimensions are smaller than this bitmap's dimensions, you can
    /// specifie
    /// an offset by which src bitmap is shifted.
    void mergeBitmap(const Bitmap<T> &src, int x_offset=0, int y_offset=0){
       for (int x=0; x<src.width_; x++){
            for (int y=0; y<src.height_; y++){
                int xdst = x+x_offset;
                int ydst = y+y_offset;
                for (int ch=0; ch<channels_; ch++){
                    setElem(xdst, ydst,
                        src.data_[(src.width_*y + x)*src.channels_ + ch], ch);
                }
            }
       }
    }

    /// Performs a maximum blending of this bitmap and src bitmap:
    /// The last channel of the bitmaps is considered the alpha channel.
    /// For each pixel either the source or destination pixel is chosen,
    /// depending on whose alpha value is greater. If src's dimensions are
    /// smaller than this bitmap's dimensions, you can specifie
    /// an offset by which src bitmap is shifted.
    void blendMax(const Bitmap<T> &src, int x_offset=0, int y_offset=0){
        for (int x=0; x<src.width_; x++){
            for (int y=0; y<src.height_; y++){
                int xdst = x+x_offset;
                int ydst = y+y_offset;
                if (getElem(xdst, ydst, channels_-1) <
                    src.data_[(src.width_*y+x)*src.channels_+(src.channels_-1)] ){
                        for (int ch=0; ch<channels_; ch++){
                            setElem(xdst, ydst,
                                src.data_[(src.width_*y + x)*src.channels_ + ch],
                                ch);
                        }
                }
            }
        }
    }

    /// Performs an alpha blending of this bitmap and the source bitmap.
    /// The last channel is considered the alpha channel.
    ///     source_alpha = source alpha
    ///     dest_alpha = (1-source_alpha)*dest_alpha
    /// If src's dimensions are smaller than this bitmap's dimensions, you can specifie
    /// an offset by which src bitmap is shifted.
    void blend(const Bitmap<T> &src, int x_offset=0, int y_offset=0){
        for (int x=0; x<src.width_; x++){
            for (int y=0; y<src.height_; y++){
                int xdst = x+x_offset;
                int ydst = y+y_offset;
                T srcAlpha = src.data_[(src.width_*y+x)*src.channels_+(src.channels_-1)];
                T dstAlpha = (1-srcAlpha)*getElem(xdst, ydst, 3);
                for (int ch=0; ch<channels_; ch++){
                    T elem = getElem(xdst, ydst, ch)*dstAlpha +
                        src.data_[(src.width_*y + x)*src.channels_ + ch]*srcAlpha;
                    setElem(xdst, ydst, elem, ch);
                }
                setElem(xdst, ydst, srcAlpha + (1-srcAlpha)*dstAlpha, channels_-1);
            }
        }
    }

    /// Deletes the bitmap's data array.
    void freeData(){
        delete[] data_;
        data_ = 0;
    }

    /// Returns the bitmap's width
    inline int getWidth(){
        return width_;
    }
    /// Returns the bitmap's height
    inline int getHeight(){
        return height_;
    }

    /// Sets image's dimension to width, height
    /// image data is cleared!
    void setDim(int width, int height, int channels=1){
        delete[] data_;
        data_ = new T[width*height*channels];
        width_ = width;
        height_ = height;
        channels_ = channels;
    }

protected:
    T* data_;
    int width_;
    int height_;
    int channels_;
};

//----------------------------------------------------------------------------

/// A specialization of Bitmap for first hit positions.
/// Modifies / adds some functions useful for first hit positions.
class FirstHitBitmap : public Bitmap<float> {

public:

    /// Returns the first hit point at (x,y) in object! coordinates.
    tgt::vec3 getPixel3Ch(int x, int y){
        return Bitmap<float>::getPixel3Ch(x, y)*2.f*volumeSize_ - volumeSize_;
    }

    /// Calculates the surface normal at (x,y) by using its neighbouring
    /// first hit points. The normal is directed to the viewer and
    /// calculated in object coordinates.
    tgt::vec3 calcNormal(int x, int y){
        tgt::vec3 xtangent = getPixel3Ch(x+1,y) - getPixel3Ch(x-1,y);
        tgt::vec3 ytangent = getPixel3Ch(x,y+1) - getPixel3Ch(x,y-1);
        return tgt::normalize( tgt::cross(xtangent, ytangent) );
    }

    /// Transforms a point from object coordinates of the proxy geometry to
    /// viewport coordinates.
    /// @param round if set to true, the x and y components are rounded
    /// to the next integer
    tgt::vec3 projectToViewport(tgt::vec3 point, bool round=true){
        tgt::vec4 p4D = transformationMatrix_*tgt::vec4(point, 1.f);
        p4D /= p4D.w; // perspective correction
        tgt::vec3 result;
        if (round)
            result = tgt::vec3( tgt::round(p4D.x), tgt::round(p4D.y) , p4D.z);
        else
            result = tgt::vec3( p4D.x, p4D.y, p4D.z );

        return result;
    }

    /// Transforms a point from viewport coordinates to object coordinates.
    tgt::vec3 inverseProject(tgt::vec3 point){
        tgt::vec4 p4D = transformationMatrixInverse_*tgt::vec4(point, 1.f);
        p4D /= p4D.w; // perspective correction
        tgt::vec3 result = tgt::vec3(p4D.x, p4D.y, p4D.z);
        return result;
    }

    /// Sets the dimension of the proxy geometry.
    void setVolumeSize(tgt::vec3 volumeSize){
        volumeSize_ = volumeSize;
    }

    /// Calculates the matrix that transforms from object coordinates to viewport
    /// coordinates and saves it internally.
    void calcTransformationMatrix(tgt::mat4 viewMatrix,
            tgt::mat4 projectionMatrix,
            tgt::ivec2 viewportSize) {
        tgt::mat4 scalemat = tgt::mat4::createScale(
                tgt::vec3(viewportSize.x/2.f, viewportSize.y/2.f, 1.f/2.f) );
        tgt::mat4 transmat = tgt::mat4::createTranslation( tgt::vec3(1.f) );

        transformationMatrix_ = scalemat*transmat*projectionMatrix*viewMatrix;
        transformationMatrix_.invert(transformationMatrixInverse_);
    }

private:

    /// Scale factors for scaling from texture coordinates
    /// to the proxy geometry's object coordinates.
    tgt::vec3 volumeSize_;

    /// Matrices for transforming object coordinates to viewport coordinates
    /// and vice versa.
    tgt::mat4 transformationMatrix_;
    tgt::mat4 transformationMatrixInverse_;

};

} // namespace labeling

#endif //VRN_LABELINGMATH_H
