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

#ifndef VRN_REALWORLDMAPPINGMETADATA_H
#define VRN_REALWORLDMAPPINGMETADATA_H

#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/datastructures/volume/volumeelement.h"
#include "templatemetadata.h"

namespace voreen {

/// Stores information to map from normalized values (@see VolumeRAM::getVoxelNormalized) to real world (e.g., HU) values.
class VRN_CORE_API RealWorldMapping : public Serializable {
public:
    RealWorldMapping();
    RealWorldMapping(float scale, float offset, std::string unit);
    RealWorldMapping(tgt::vec2 range, std::string unit);

    /// Testing for equality with another Real World Mapping. Evaluates to true if offset, scale and unit are equal.
    bool operator==(const RealWorldMapping& other) const;
    /// Testing for inequality with another Real World Mapping. Evaluates to true if any of offset, scale and unit differs.
    bool operator!=(const RealWorldMapping& other) const;

    virtual void serialize(Serializer& s) const; ///< @see Serializable::serialize
    virtual void deserialize(Deserializer& s);   ///< @see Serializable::deserialize

    /**
     * @brief Transforms normalized value to real world
     * @return (normalized * scale) + offset
     */
    float normalizedToRealWorld(float normalized) const;
    /**
     * @brief Transforms real world value to normalized
     * @return (realWorld - offset) / scale
     */
    float realWorldToNormalized(float realWorld) const;

    /**
     * Returns the inverted mapping.
     */
    RealWorldMapping getInverseMapping() const;

    std::string getUnit() const;
    void setUnit(std::string unit);
    float getScale() const;
    void setScale(float scale);
    float getOffset() const;
    void setOffset(float offset);

    /// Create mapping from normalized to original values.
    template<typename T>
    static RealWorldMapping createDenormalizingMapping() {
        float scale = 1.0f;
        float offset = 0.0f;

        if(VolumeElement<T>::isInteger()) {
            if(VolumeElement<T>::isSigned()) {
                scale = (VolumeElement<T>::rangeMaxElement() - VolumeElement<T>::rangeMinElement()) / 2.0f;
                offset = (VolumeElement<T>::rangeMaxElement() + VolumeElement<T>::rangeMinElement()) / 2.0f;
            }
            else {
                scale = VolumeElement<T>::rangeMaxElement();
            }
        }

        return RealWorldMapping(scale, offset, "");
    }

    /**
     * Create mapping from normalized to original values.
     * Can be used instead of the templated version if the format is not known at compile time.
     * @param format Voreen style string base type identifier
     * @throws tgt::UnsupportedFormatException in case an unknown format is given as parameter
     */
    static RealWorldMapping createDenormalizingMapping(const std::string& baseType);

    /// Create a combined mapping, first applying m1, then m2 (using unit from m2).
    static RealWorldMapping combine(const RealWorldMapping& m1, const RealWorldMapping& m2) {
        return RealWorldMapping(m1.getScale()*m2.getScale(), (m2.getScale()*m1.getOffset()) + m2.getOffset(), m2.getUnit());
    }

private:
    float scale_;
    float offset_;

    std::string unit_;
};

inline std::ostream& operator<<(std::ostream& s, const RealWorldMapping& rwm) {
    return (s << "Scale: " << rwm.getScale() << " Offset: " << rwm.getOffset() << " Unit: " << rwm.getUnit());
}

///Metadata encapsulating RealWorldMapping
class VRN_CORE_API RealWorldMappingMetaData : public TemplateMetaData<RealWorldMapping> {
public:
    RealWorldMappingMetaData() : TemplateMetaData<RealWorldMapping>() {}
    RealWorldMappingMetaData(RealWorldMapping value) : TemplateMetaData<RealWorldMapping>(value) {}
    RealWorldMappingMetaData(float slope, float offset, std::string unit);

    virtual MetaDataBase* clone() const      { return new RealWorldMappingMetaData(getValue()); }
    virtual std::string getClassName() const { return "RealWorldMappingMetaData"; }
    virtual MetaDataBase* create() const { return new RealWorldMappingMetaData(); }

    virtual std::string toString() const;
    virtual std::string toString(const std::string& component) const;
};

} // namespace

#endif
