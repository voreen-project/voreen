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

#ifndef VRN_STREAMLINE_H
#define VRN_STREAMLINE_H

#include "voreen/core/voreencoreapi.h"

#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/io/serialization/xmlserializer.h"
#include "voreen/core/io/serialization/xmldeserializer.h"

#include "tgt/vector.h"

#include <deque>

namespace voreen {

    /**
     * Datastructure used to represent streamlines. It is used in the flowreen module.
     * A streamline consists of multiple StreamlineElements each storing a position and the velocity at this position.
     */
class VRN_CORE_API Streamline : public Serializable {
    friend class StreamlineList;
public:
    /**
     * A streamline consits of a vector of these elements.
     */
    struct StreamlineElement {
        tgt::vec3 position_;    ///< position in voxel space
        tgt::vec3 velocity_;

        StreamlineElement()
            : position_(tgt::vec3::zero)
            , velocity_(tgt::vec3::zero)
        {}

        StreamlineElement(tgt::vec3 position, tgt::vec3 velocity)
            : position_(position)
            , velocity_(velocity)
        {
            //tgtAssert(tgt::length(velocity) > 16.f, "ohoh"); debug
        }
    };

    /** Constructor */
    Streamline();
    /** Destructor */
    ~Streamline();

    //----------------
    //  Construction
    //----------------
    /** Adds and copies an element. */
    void addElementAtEnd(const StreamlineElement& element);
    /** Adds and copies an element. */
    void addElementAtFront(const StreamlineElement& element);

    //----------------
    //  Access
    //----------------
    const StreamlineElement& getElementAt(size_t pos) const;
    const StreamlineElement& getFirstElement() const;
    const StreamlineElement& getLastElement() const;
    size_t getNumElements() const;
    float getMinMagnitude() const;
    float getMaxMagnitude() const;
    float getLength() const;

    //----------------
    //  Utility
    //----------------
    /** Resamples this Streamline to a similar one consisting of the specified amount of elements. */
    Streamline resample(size_t samples) const;

    //----------------
    //  Storage
    //----------------
    /** Used to save as CSV file. Transforms the voxel position to world space. */
    std::string toCSVString(const tgt::mat4& transfomationMatrix = tgt::mat4::identity,
                            const tgt::mat4& velocityTransfomationMatrix = tgt::mat4::identity) const;
    /** @override */
    virtual void serialize(Serializer& s) const;
    /** @override */
    virtual void deserialize(Deserializer& s);

    //----------------
    //  Members
    //----------------
private:
   std::deque<StreamlineElement> streamlineElements_;   ///< list of all streamline elements from front to back
   float minMagnitude_;                                 ///< min magnitude value for color map configuration
   float maxMagnitude_;                                 ///< max magnitude value for color map configuration
   float length_;                                       ///< the actual length of the streamline
};

}   // namespace

#endif
