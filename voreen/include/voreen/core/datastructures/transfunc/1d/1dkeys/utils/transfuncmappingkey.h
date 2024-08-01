/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_TRANSFUNCMAPPINGKEY_H
#define VRN_TRANSFUNCMAPPINGKEY_H

#include "voreen/core/voreenobject.h"
#include "voreen/core/io/serialization/serialization.h"

#include "tgt/vector.h"
#include "tgt/texture.h"

namespace voreen {

/**
 * One of multiple keys that define a TransfuncIntensity.
 * Each key has an intensity at which it is located and a color.
 * Furthermore it can be split in two parts with two different colors.
 */
class VRN_CORE_API TransFuncMappingKey : public VoreenSerializableObject {
public:
    /**
     * Constructor
     *
     * @param i intensity at which this key is located
     * @param color the color of this key
     */
    TransFuncMappingKey(float i, const tgt::col4& color);

    /**
     * Standard Destructor
     */
    ~TransFuncMappingKey();

    virtual std::string getClassName() const    { return "TransFuncMappingKey";     }
    virtual TransFuncMappingKey* create() const { return new TransFuncMappingKey(); }

    /**
     * Operator to compare two TransFuncMappingKeys. True is returned when both keys
     * have the same intensity, color and split status. Otherwise false is returned.
     *
     * @param key key that is compared with this key
     * @return true if both keys have same split status, same color and same intensity, false otherwise
     */
    bool operator==(const TransFuncMappingKey& key);

    /**
     * Operator to compare two TransFuncMappingKeys. True is returned when both keys
     * differ in color, intensity or split status. Otherwise false is returned.
     *
     * @param key key that is compared with this key
     * @return true if both keys differ in split status, color or intensity, false otherwise
     */
    bool operator!=(const TransFuncMappingKey& key);

    /**
     * Sets the color of the left part of the key to the given value.
     *
     * @param color color the left part of the key will be set to
     */
    void setColorL(const tgt::col4& color);

    /**
     * Sets the color of the left part of the key to the given value.
     *
     * @param color color the left part of the key will be set to
     */
    void setColorL(const tgt::ivec4& color);

    /**
     * Returns the color of the left part of the key.
     *
     * @return color of the left part of the key
     */
    const tgt::col4& getColorL() const;
    tgt::col4& getColorL();

    /**
     * Sets the color of the right part of the key to the given value.
     *
     * @param color color the right part of the key will be set to
     */
    void setColorR(const tgt::col4& color);

    /**
     * Sets the color of the right part of the key to the given value.
     *
     * @param color color the right part of the key will be set to
     */
    void setColorR(const tgt::ivec4& color);

    /**
     * Returns the color of the right part of the key.
     *
     * @return color of the right part of the key
     */
    const tgt::col4& getColorR() const;
    tgt::col4& getColorR();

    /**
     * Sets the alpha value of the right part of the key to the given value.
     * The value has to be in the range 0.f-1.f.
     *
     * @param a alpha value the right part of the key will be set to
     */
    void setAlphaR(float a);

    /**
     * Sets the alpha value of the left part of the key to the given value.
     * The value has to be in the range 0.f-1.f.
     *
     * @param a alpha value the left part of the key will be set to
     */
    void setAlphaL(float a);

    /**
     * Returns the alpha value of the right part of the key as float (0.f-1.f).
     *
     * @return alpha value of the right part of the key
     */
    float getAlphaR();

    /**
     * Returns the alpha value of the left part of the key as float (0.f-1.f).
     *
     * @return alpha value of the left part of the key
     */
    float getAlphaL();

    /**
     * Returns whether the key is split or not.
     *
     * @return true if the key is split, false otherwise.
     */
    bool isSplit() const;

    /**
     * Splits or unsplits this key.
     *
     * @param split true to split, false to merge
     * @param useLeft in case of joining: use the left color (else use right)
     */
    void setSplit(bool split, bool useLeft=true);

    /**
     * Returns the intensity at which the key is located.
     *
     * @return intensity at which the key is located
     */
    float getIntensity() const;

    /**
     * Sets the intensity of the key to the given value. The keys of the transfer function
     * have to be sorted after this.
     *
     * @param i new intensity of the key
     */
    void setIntensity(float i);

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
    TransFuncMappingKey* clone() const;

private:
    friend class XmlDeserializer;

    /**
     * Default constructor needed for serialization purposes.
     */
    TransFuncMappingKey();

    float intensity_;  ///< intensity at which the key is located
    tgt::col4 colorL_; ///< color of the left part of the key
    tgt::col4 colorR_; ///< color of the right part of the key
    bool split_;       ///< is the key split?
};

} // namespace voreen

#endif // VRN_TRANSFUNCMAPPINGKEY_H
