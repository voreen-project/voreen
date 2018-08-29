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

#ifndef VRN_COLORMAP_H
#define VRN_COLORMAP_H

#include "voreen/core/voreencoreapi.h"
#include "voreen/core/io/serialization/serializable.h"
#include "tgt/vector.h"

#include <vector>
#include <iterator>

namespace voreen {

/**
 * \brief   This class encapsulates \sa tgt::Color objects and provides methods for automatic color creation.
 *
 * A ColorMap object can only be instantiated via the static create methods and holds at least two colors.
 * Colors of this ColorMap can be aquired via local getter function or Iterators. InterpolationIterator will
 * interpolate between the held colors in n steps (n given at iterator instantiation). GeneratingIterator
 * will create a infinite sequence of colors generated out of the Colors in this ColorMap.
 **/
class VRN_CORE_API ColorMap : public Serializable {

// BoxObject needs a parameterless constructor to work with ColorMap, though we grant it access to our
// private default constructor
friend class BoxObject;

public:
    /**
     * Interpolates between the held colors of a given ColorMap in n steps (color map and n given
     * at iterator instantiation).
     **/
    class InterpolationIterator : public std::iterator<std::forward_iterator_tag, tgt::Color> {
    public:
        /// constructor
        InterpolationIterator(const ColorMap* map, int currentIndex, int maximum_)
        : map_(map)
        , currentIndex_(currentIndex)
        , maximum_(maximum_)
        {}

        /// copy constructor
        InterpolationIterator(const InterpolationIterator& rhs)
        : map_(rhs.map_)
        , currentIndex_(rhs.currentIndex_)
        , maximum_(rhs.maximum_)
        {}

        InterpolationIterator& operator++() {
            ++currentIndex_; return *this;
        }

        bool operator==(const InterpolationIterator& rhs) {
            return currentIndex_ == rhs.currentIndex_;
        }
        bool operator!=(const InterpolationIterator& rhs) {
            return currentIndex_ != rhs.currentIndex_;
        }

        tgt::Color operator*() {
            return map_->getColorAtPosition(static_cast<float>(currentIndex_)/maximum_);
        }

    private:
        const ColorMap* map_;

        int currentIndex_;
        int maximum_;
    };

    /**
     * Creates an infinite sequence of colors generated out of the colors in the given ColorMap.
     **/
    class GeneratingIterator : public std::iterator<std::forward_iterator_tag, tgt::Color> {
        public:
            /// constructor
            GeneratingIterator(const ColorMap* map, int currentIndex)
            : map_(map)
            , currentIndex_(currentIndex)
            { }

            /// copy constructor
            GeneratingIterator(const GeneratingIterator& rhs)
            : map_(rhs.map_)
            , currentIndex_(rhs.currentIndex_)
            {}

            GeneratingIterator& operator++() {
                ++currentIndex_; return *this;
            }

            bool operator==(const GeneratingIterator& rhs) {
                return currentIndex_ == rhs.currentIndex_;
            }
            bool operator!=(const GeneratingIterator& rhs) {
                return currentIndex_ != rhs.currentIndex_;
            }


            tgt::Color operator*();
        private:
            const ColorMap* map_;

            int currentIndex_;
    };

    friend class InterpolationIterator;
    friend class GeneratingIterator;

    /**
     * Constructor is private as ColorMaps shall only be instantiated via static
     * create Methods
     **/
    ColorMap();

    // predefined color maps
    // please change getColorMapLabels and createColorMap if you add a colormap
    static ColorMap createSpectral();
    static ColorMap createColdHot();
    static ColorMap createTango();
    static ColorMap createVisifire1();
    static ColorMap createVisifire2();
    static ColorMap createCandleLight();
    static ColorMap createSandyShades();
    static ColorMap createGermany();

    // returns names of all predefined colormaps
    static std::vector<std::string> getColorMapLabels();

    //creates the predefinied colormap with index (in getColorMapLabels)
    static ColorMap createColorMap(int index);
    static ColorMap createSequential(const tgt::Color& base);
    static ColorMap createFromVector(const std::vector<tgt::Color>& colors);

    /**
    * @see Property::serialize
    */
    virtual void serialize(Serializer& s) const;

    /**
    * @see Property::deserialize
    */
    virtual void deserialize(Deserializer& s);

    std::string toString() const;

    void setName(const std::string& name);

    /**
     * \brief adds the color \a color to the end of this ColorMap
     *
     * \param   color   color to add to the ColorMap
     **/
    void addColorLast(const tgt::Color& color);

    /**
     * \brief adds the color \a color at position \a position into this ColorMap
     *
     * \param   color       color to add to the ColorMap
     * \param   position    0-based index where to add the color into the ColorMap
     **/
    void addColorAtIndex(const tgt::Color& color, int position);

    /**
     * \brief returns the color with index \a position
     *
     * \param   position    0-based index which color to return
     **/
    tgt::Color getColorAtIndex(int position) const;

    /**
     * \brief   returns the color at position \a position within interval [0,1].
     * If no such color is defined ColorMap will do linear interpolation between surrounding colors.
     *
     * \param   position    position of which color to return, must be within interval [0,1]
     **/
    tgt::Color getColorAtPosition(float position) const;

    int getColorCount() const;

    InterpolationIterator getInterpolationIterator(int count) const;

    GeneratingIterator getGeneratingIterator() const;

    bool operator!=(const ColorMap& rhs) const;
    bool operator==(const ColorMap& rhs) const;

private:
    static tgt::vec4 rgbToHsv(tgt::Color rgb);

    static tgt::Color hsvToRgb(tgt::vec4 hsv);

    std::vector<tgt::Color> colors_;
    std::string name_;

};

} // namespace voreen

#endif // VRN_COLORMAP_H

