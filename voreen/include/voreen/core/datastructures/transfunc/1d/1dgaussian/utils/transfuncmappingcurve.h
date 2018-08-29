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

#ifndef VRN_TRANSFUNCMAPPINGCURVE_H
#define VRN_TRANSFUNCMAPPINGCURVE_H

#include "voreen/core/voreenobject.h"
#include "voreen/core/io/serialization/serialization.h"

#include "tgt/vector.h"
#include "tgt/texture.h"

namespace voreen {

/*
* Enum to identify the separate parts of a curve.
* Mainly used in TransFuncMappingCurve::setUnicolor()
*/
enum CurvePart {
    BASE_LEFT = 0,
    BASE_RIGHT = 1,
    PEAK_LEFT = 2,
    PEAK_RIGHT = 3,
    NO_KEY = -1
};

/**
 * One of the curves that define a TransfuncIntensity.
 * Each curve has an intensity at which it is located and a color.
 * Furthermore it can be split in two parts with two different colors.
 */
class VRN_CORE_API TransFuncMappingCurve : public VoreenSerializableObject {
public:
    /**
     * Constructor
     *
     * @param i intensity at which the peak of this curve is located
     * @param width the variance of this curve's Gauss function
     * @param baseValue the base value of this curve's Gauss function
     * @param color the color of this curve
     */
    TransFuncMappingCurve(float i, float width, float baseValue, const tgt::col4& color);

    /**
     * Standard Destructor
     */
    ~TransFuncMappingCurve();

    virtual std::string getClassName() const    { return "TransFuncMappingCurve";     }
    virtual TransFuncMappingCurve* create() const { return new TransFuncMappingCurve(); }

    /**
     * Operator to compare two TransFuncMappingCurves. True is returned when both curves
     * have the same intensity, width, base value, color and split status. Otherwise false is returned.
     *
     * @param curve curve that is compared with this curve
     * @return true if both curves have same attributes, false otherwise
     */
    bool operator==(const TransFuncMappingCurve& curve);

    /**
     * Operator to compare two TransFuncMappingCurves. True is returned when both curves
     * differ in color, intensity, width, base value or split status. Otherwise false is returned.
     *
     * @param curve curve that is compared with this curve
     * @return true if both curves have the same attributes, false otherwise
     */
    bool operator!=(const TransFuncMappingCurve& curve);

    /**
     * Sets the color of the left part of the peak key to the given value.
     *
     * @param color color the left part of the peak key will be set to
     */
    void setColorL(const tgt::col4& color);

    /**
     * Sets the color of the left part of the peak key to the given value.
     *
     * @param color color the left part of the peak key will be set to
     */
    void setColorL(const tgt::ivec4& color);

    /**
     * Returns the color of the left part of the peak key.
     *
     * @return color of the left part of the peak key
     */
    const tgt::col4& getColorL() const;
    tgt::col4& getColorL();

    /**
     * Sets the color of the right part of the peak key to the given value.
     *
     * @param color color the right part of the peak key will be set to
     */
    void setColorR(const tgt::col4& color);

    /**
     * Sets the color of the right part of the peak key to the given value.
     *
     * @param color color the right part of the peak key will be set to
     */
    void setColorR(const tgt::ivec4& color);

    /**
     * Returns the color of the right part of the peak key.
     *
     * @return color of the right part of the peak key
     */
    const tgt::col4& getColorR() const;
    tgt::col4& getColorR();

    /**
    * Sets the color of the left base key to the given value.
    *
    * @param color color the left base key will be set to
    */
    void setBaseColorL(const tgt::col4& color);

    /**
    * Sets the color of the left base key to the given value.
    *
    * @param color color the left base key will be set to
    */
    void setBaseColorL(const tgt::ivec4& color);

    /**
    * Returns the color of the left base key.
    *
    * @return color of the left base key
    */
    const tgt::col4& getBaseColorL() const;
    tgt::col4& getBaseColorL();

    /**
    * Sets the color of the right base key to the given value.
    *
    * @param color color the right base key will be set to
    */
    void setBaseColorR(const tgt::col4& color);

    /**
    * Sets the color of the right base key to the given value.
    *
    * @param color color the right base key will be set to
    */
    void setBaseColorR(const tgt::ivec4& color);

    /**
    * Returns the color of the right base key.
    *
    * @return color of the right base key
    */
    const tgt::col4& getBaseColorR() const;
    tgt::col4& getBaseColorR();

    /**
     * Sets the alpha value of the right part of the peak key to the given value.
     * The value has to be in the range 0.f-1.f.
     *
     * @param a alpha value the right part of the peak key will be set to
     */
    void setAlphaR(float a);

    /**
     * Sets the alpha value of the left part of the peak key to the given value.
     * The value has to be in the range 0.f-1.f.
     *
     * @param a alpha value the left part of the peak key will be set to
     */
    void setAlphaL(float a);

    /**
     * Returns the alpha value of the right part of the peak key as float (0.f-1.f).
     *
     * @return alpha value of the right part of the peak key
     */
    float getAlphaR();

    /**
     * Returns the alpha value of the left part of the peak key as float (0.f-1.f).
     *
     * @return alpha value of the left part of the peak key
     */
    float getAlphaL();

    /**
     * Returns the intensity at which the curve is located.
     *
     * @return intensity at which the curve is located
     */
    float getIntensity() const;

    /**
     * Sets the intensity of the curve to the given value.
     *
     * @param i new intensity of the key
     */
    void setIntensity(float i);

    /**
     * Returns the variance of this Gauss function.
     *
     * @return width of this curve's Gauss function
     */
    float getWidth() const;

    /**
     * Sets the variance of this Gauss function to the given value.
     *
     * @param w new width of this curve's Gauss function
     */
    void setWidth(float w);

    /**
    * Returns the opacity value to which this curves's Gauss function converges at the borders.
    * The peak of the function can be below as well as above the base value.
    *
    * @return base opacity of the curve
    */
    float getBaseValue() const;

    /**
    * Sets the opacity value to which this curve's Gauss function converges at the borders.
    * The peak of the function can be below as well as above the base value.
    *
    * @param w new base opacity of the curve
    */
    void setBaseValue(float a);

    /**
    * Returns the opacity value at the given intensity resulting from this curve.
    * The return value lies between 0 and 1 (inclusive).
    *
    * Based on the function: f(x) = e^(-(x - intensity)²/(2 * width)) * (alpha - baseValue) + baseValue
    *
    * @param opacity at intensity (between 0 and 1)
    */
    float getOpacityAt(float intensity);

    /**
    * Returns the color value with the appropriate alpha value at intensity
    * resulting from this curves Gauss function.
    * Note that especially the alpha value of the returned color will be
    * much more imprecise than the float value returned by getOpacityAt().
    *
    * @return color value at intensity with rgba in [0, 255]
    */
    tgt::col4 getColorAt(float intensity);

    /**
    * Returns whether the curve is split or not.
    *
    * @return true if the curve is split, false otherwise.
    */
    bool isSplit() const;

    /**
    * Splits or unsplits this curve.
    *
    * @param split true to split, false to merge
    * @param useLeft in case of joining: use the left color (else use right)
    */
    void setSplit(bool split, bool useLeft = true);

    /**
    * Sets the curve to only use one color for all three keys (peak, left/right base)
    * or not. If different key colors were set before, the color of the specified key
    * determines the new overall color. In an unicolor curve all setter methods
    * for the colors change all three key colors (rgb) equally,
    * but without changing their alpha values.
    *
    * @param unicolor if the same color is used for alle three keys of this curve
    * @param key whose color is transferred to all other keys when the curve is set to unicolor
    */
    void setUnicolor(bool unicolor, CurvePart relevantKey = PEAK_LEFT);

    /**
    * Returns if the curve uses only one color for all three keys (peak, left/right base)
    * or not. In an unicolor curve all setter methods for the colors change all three
    * key colors (rgb) equally, but without changing their alpha values.
    *
    * @return if the same color is used for alle three keys of this curve
    */
    bool isUnicolor();

    /**
    * Activates or deactivates this curve.
    * Only active curves shall be used in the calculation of the transfer function.
    *
    * @param active new state of this curve
    */
    void setActive(bool active);

    /**
    * Returns wether this curve is activated or not.
    * Only active curves shall be used in the calculation of the transfer function.
    *
    * @return true when this curve is active
    */
    bool isActive();

    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);

    /**
     * Returns a copy of this object.
     */
    TransFuncMappingCurve* clone() const;

private:
    friend class XmlDeserializer;

    /**
    * Sets all base and peak colors to the given value.
    * All alpha values remain unchanged if setAlpha is false.
    *
    * @param color the target color for all color attributes
    * @param setAlpha if the alpha values will also be set
    */
    void setAllColors(const tgt::col4& color, bool setAlpha = false);

    /**
     * Default constructor needed for serialization purposes.
     */
    TransFuncMappingCurve();

    float intensity_;  ///< intensity at which the curve is located
    float width_;      ///< the variance of the Gauss function represented by this curve
    float baseValue_;  ///< the basse opaicty value of this curve
    tgt::col4 colorL_; ///< color of the left part at the peak of the curve.
                       ///  The alpha value of this color determines the curves base value.
    tgt::col4 colorR_; ///< color of the right part at the peak of the curve
    tgt::col4 baseColorL_;    ///< color of the left part at the base line of the curve
    tgt::col4 baseColorR_;    ///< color of the right part at the base line of the curve
    bool split_;       ///< is the curve split?
    bool unicolor_;    ///< is the color for all three keys (peak, left/right bottom) the same?
    bool active_;      ///< is the curve active?
};

} // namespace voreen

#endif // VRN_TRANSFUNCMAPPINGCURVE_H
