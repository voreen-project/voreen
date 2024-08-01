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

#ifndef VRN_STREAMLINE_H
#define VRN_STREAMLINE_H

#include "voreen/core/voreencoreapi.h"

#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/io/serialization/xmlserializer.h"
#include "voreen/core/io/serialization/xmldeserializer.h"

#include "voreen/core/utils/statistics.h"

#include "tgt/vector.h"

#include <deque>

namespace voreen {

/**
 * Datastructure used to represent streamlines. It is used in the flowanalysis module.
 * A streamline consists of multiple StreamlineElements each storing a position and the velocity at this position.
 */
class VRN_CORE_API Streamline : public Serializable {
public:

    /**
     * A streamline consists of a vector of these elements.
     */
    struct StreamlineElement {
        tgt::vec3 position_;    ///< position in world space
        tgt::vec3 velocity_;    ///< local velocity in mm/s
        float radius_;          ///< local radius in world space, primarily used for bundles
        float time_;            ///< time point in seconds, primarily used for pathlines

        StreamlineElement()
            : position_(tgt::vec3::zero)
            , velocity_(tgt::vec3::zero)
            , radius_(0.0f)
            , time_(0.0f)
        {}

        StreamlineElement(const tgt::vec3& position, const tgt::vec3& velocity, float radius = 0.0f, float time = 0.0f)
            : position_(position)
            , velocity_(velocity)
            , radius_(radius)
            , time_(time)
        {}
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

    //----------------
    //  Utility
    //----------------
    /** Resamples this Streamline to a similar one consisting of the specified amount of elements. */
    Streamline resample(size_t samples) const;

    /** Returns statistics of the angle (in radians) between to consecutive elements. */
    const Statistics& getCurvatureStatistics() const;

    /** Returns statistics of the magnitude of all elements. */
    const Statistics& getMagnitudeStatistics() const;

    /** Returns the minimum magnitude of all elements. */
    float getMinMagnitude() const;

    /** Returns the maximum magnitude of all elements. */
    float getMaxMagnitude() const;

    /** Returns the length in physical space. */
    float getPhysicalLength() const;

    /** Returns the temporal range covered by the streamline. */
    tgt::vec2 getTemporalRange() const;

    //----------------
    //  Storage
    //----------------
    /** Used to save as CSV file. */
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
    Statistics magnitudeStatistics_;                     ///< statistics of the contained magnitudes
    Statistics curvatureStatistics_;                     ///< statistics of the lines curvature
    float physicalLength_;                               ///< total physical length (in mm)
};

}   // namespace

#endif
