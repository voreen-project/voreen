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

#ifndef VRN_ANIMATEDPROCESSOR_H
#define VRN_ANIMATEDPROCESSOR_H

#include <string>
#include <vector>
#include "voreen/core/animation/propertytimeline.h"
#include "voreen/core/io/serialization/serialization.h"

namespace voreen {

class Processor;
class Property;

/**
 * This class manages all animations for one processor.
 */
class VRN_CORE_API AnimatedProcessor : public Serializable {

public:
    /**
     * Constructor, gets a pointer to the processor which should be animated.
     */
    AnimatedProcessor(Processor* processor);
    ~AnimatedProcessor();

    /**
     * Returns an array of all propertytimelines which belong to this processor.
     */
    const std::vector<PropertyTimeline*>& getPropertyTimelines() const;

    /**
     * Renders the complete Animation at a given time - respectively sets all properties of all processors to the value defined in the animation
     */
    void renderAt(float time);

    /**
     * Returns a pointer to the corresponding processor
     */
    const Processor* getCorrespondingProcessor() const;
    Processor* getCorrespondingProcessor();

    /**
     * Returns the name of the animated processor.
     */
    const std::string getProcessorName() const;

    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);

    bool isPropertyAnimated(const Property* p);

    /**
     * Adds a new timeline for the given property.
     */
    PropertyTimeline* addTimeline(Property* prop);

    /**
     * Removes the property timeline from the animated processor and returns a pointer to it for deleting it later on.
     */
    void removeTimeline(Property* prop);

protected:
    static const std::string loggerCat_; ///< category used in logging

private:
    friend class XmlDeserializer;

    AnimatedProcessor();

    Processor* processor_;
    std::vector<PropertyTimeline*> properties_;
};

} // namespace voreen

#endif
