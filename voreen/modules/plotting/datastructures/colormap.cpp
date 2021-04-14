/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "colormap.h"
#include "voreen/core/io/serialization/serialization.h"

#include <cmath>

namespace voreen {

tgt::Color ColorMap::GeneratingIterator::operator*() {
    int cc = map_->getColorCount();
    tgt::Color base = map_->getColorAtIndex(currentIndex_ % cc);

    if (currentIndex_ > 0 && currentIndex_ >= cc) {
        // to create new colors we turn the hue by 37 degrees and randomly modify saturation and light
        tgt::vec4 foo = ColorMap::rgbToHsv(base);

        float mult = static_cast<float>(floor(static_cast<float>(currentIndex_)/cc));
        foo.x += mult*37.f;
        foo.x = fmod(foo.x, 360.f);

        float mult2 = (base.r + base.g + base.b)/6.f;
        foo.y *= (mult2 + 1.f);

        // normalize back to [0.3, 0.8]
        while (foo.y > .8f)
            foo.y -= .5f;

        base = ColorMap::hsvToRgb(foo);
    }

    return base;
}

ColorMap ColorMap::createSpectral() {
    ColorMap cm;
    cm.colors_.clear();
    cm.addColorLast(tgt::Color(1.f, 0.f, 0.f, 1.f));
    cm.addColorLast(tgt::Color(1.f, 1.f, 0.f, 1.f));
    cm.addColorLast(tgt::Color(0.f, 1.f, 0.f, 1.f));
    cm.addColorLast(tgt::Color(0.f, 0.f, 1.f, 1.f));
    cm.setName("Spectral");
    return cm;
}

ColorMap ColorMap::createColdHot() {
    ColorMap cm;
    cm.colors_.clear();
    cm.addColorLast(tgt::Color(0.f, 0.f, 1.f, 1.f));
    cm.addColorLast(tgt::Color(1.f, 1.f, 1.f, 1.f));
    cm.addColorLast(tgt::Color(1.f, 0.f, 0.f, 1.f));
    cm.setName("Cold Hot");
    return cm;
}

ColorMap ColorMap::createTango() {
    ColorMap cm;
    cm.colors_.clear();
    cm.addColorLast(tgt::Color(.7656, .625, 0, 1));      //Butter
    cm.addColorLast(tgt::Color(.6172, .6016, .0234, 1)); // Chameleon
    cm.addColorLast(tgt::Color( .125, .2891, .5273, 1)); // Sky Blue
    cm.addColorLast(tgt::Color(.3594,  .207, .3984, 1)); // Plum
    cm.addColorLast(tgt::Color(.6406, 0, 0, 1));         // Scarlet Red
    cm.addColorLast(tgt::Color(.5586, .3477, .0078, 1)); // Chocolate
    cm.setName("Tango");
    return cm;
}

ColorMap ColorMap::createVisifire1() {
    ColorMap cm;
    cm.colors_.clear();
    cm.addColorLast(tgt::Color (0.26171875, 0.5234375, 0.84375, 1) );
    cm.addColorLast(tgt::Color (0.99609375, 0.6015625, 0.1796875, 1) );
    cm.addColorLast(tgt::Color (0.85546875, 0.265625, 0.24609375, 1) );
    cm.addColorLast(tgt::Color (0.65625, 0.828125, 0.30859375, 1) );
    cm.addColorLast(tgt::Color (0.51953125, 0.375, 0.69921875, 1) );
    cm.addColorLast(tgt::Color (0.234375, 0.74609375, 0.88671875, 1) );
    cm.addColorLast(tgt::Color (0.68359375, 0.84375, 0.96875, 1) );
    cm.addColorLast(tgt::Color (0.0, 0.5546875, 0.5546875, 1) );
    cm.addColorLast(tgt::Color (0.54296875, 0.7265625, 0.0, 1) );
    cm.addColorLast(tgt::Color (0.9765625, 0.73828125, 0.05859375, 1) );
    cm.addColorLast(tgt::Color (0.9765625, 0.4296875, 0.2734375, 1) );
    cm.addColorLast(tgt::Color (0.61328125, 0.03125, 0.05078125, 1) );
    cm.addColorLast(tgt::Color (0.62890625, 0.5234375, 0.7421875, 1) );
    cm.addColorLast(tgt::Color (0.796875, 0.3984375, 0.0, 1) );
    cm.addColorLast(tgt::Color (0.98828125, 0.7734375, 0.53515625, 1) );
    cm.setName("Visifire 1");
    return cm;
}

ColorMap ColorMap::createVisifire2() {
    ColorMap cm;
    cm.colors_.clear();
    cm.addColorLast(tgt::Color (0.0, 0.68359375, 0.6875, 1) );
    cm.addColorLast(tgt::Color (0.34375, 0.62890625, 0.15625, 1) );
    cm.addColorLast(tgt::Color (0.21484375, 0.3828125, 0.7265625, 1) );
    cm.addColorLast(tgt::Color (0.99609375, 0.3984375, 0.73046875, 1) );
    cm.addColorLast(tgt::Color (0.53125, 0.78125, 0.8359375, 1) );
    cm.addColorLast(tgt::Color (0.99609375, 0.796875, 0.0, 1) );
    cm.addColorLast(tgt::Color (0.18359375, 0.93359375, 0.93359375, 1) );
    cm.addColorLast(tgt::Color (0.96875, 0.25, 0.0, 1) );
    cm.addColorLast(tgt::Color (0.4765625, 0.44921875, 0.38671875, 1) );
    cm.addColorLast(tgt::Color (0.69921875, 0.7734375, 0.19140625, 1) );
    cm.addColorLast(tgt::Color (0.99609375, 0.58984375, 0.63671875, 1) );
    cm.addColorLast(tgt::Color (0.58203125, 0.42578125, 0.63671875, 1) );
    cm.addColorLast(tgt::Color (0.19140625, 0.71484375, 0.484375, 1) );
    cm.addColorLast(tgt::Color (0.6953125, 0.7421875, 0.74609375, 1) );
    cm.addColorLast(tgt::Color (0.09765625, 0.59765625, 0.86328125, 1) );
    cm.setName("Visifire 2");
    return cm;
}

ColorMap ColorMap::createCandleLight() {
    ColorMap cm;
    cm.colors_.clear();
    cm.addColorLast(tgt::Color (0.65625, 0.828125, 0.30859375, 1) );
    cm.addColorLast(tgt::Color (0.26171875, 0.5234375, 0.84375, 1) );
    cm.addColorLast(tgt::Color (0.99609375, 0.6015625, 0.1796875, 1) );
    cm.addColorLast(tgt::Color (0.51953125, 0.375, 0.69921875, 1) );
    cm.addColorLast(tgt::Color (0.234375, 0.74609375, 0.88671875, 1) );
    cm.addColorLast(tgt::Color (0.68359375, 0.84375, 0.96875, 1) );
    cm.addColorLast(tgt::Color (0.0, 0.5546875, 0.5546875, 1) );
    cm.addColorLast(tgt::Color (0.9765625, 0.73828125, 0.05859375, 1) );
    cm.addColorLast(tgt::Color (0.9765625, 0.4296875, 0.2734375, 1) );
    cm.addColorLast(tgt::Color (0.62890625, 0.5234375, 0.7421875, 1) );
    cm.addColorLast(tgt::Color (0.796875, 0.3984375, 0.0, 1) );
    cm.addColorLast(tgt::Color (0.99609375, 0.4921875, 0.0, 1) );
    cm.addColorLast(tgt::Color (0.17578125, 0.8203125, 0.46484375, 1) );
    cm.setName("Candle Light");
    return cm;
}

ColorMap ColorMap::createSandyShades() {
    ColorMap cm;
    cm.colors_.clear();
    cm.addColorLast(tgt::Color (0.6875, 0.6796875, 0.68359375, 1) );
    cm.addColorLast(tgt::Color (0.99609375, 0.7734375, 0.546875, 1) );
    cm.addColorLast(tgt::Color (0.5, 0.5, 0.0, 1) );
    cm.addColorLast(tgt::Color (0.83203125, 0.57421875, 0.3671875, 1) );
    cm.addColorLast(tgt::Color (0.38671875, 0.5, 0.5625, 1) );
    cm.addColorLast(tgt::Color (0.25, 0.35546875, 0.26953125, 1) );
    cm.addColorLast(tgt::Color (0.35546875, 0.34375, 0.32421875, 1) );
    cm.addColorLast(tgt::Color (0.55859375, 0.29296875, 0.1953125, 1) );
    cm.addColorLast(tgt::Color (0.5234375, 0.22265625, 0.00390625, 1) );
    cm.addColorLast(tgt::Color (0.42578125, 0.40234375, 0.26171875, 1) );
    cm.addColorLast(tgt::Color (0.44140625, 0.40625, 0.125, 1) );
    cm.addColorLast(tgt::Color (0.4921875, 0.375, 0.2734375, 1) );
    cm.setName("Sandy Shades");
    return cm;
}

ColorMap ColorMap::createGermany() {
    ColorMap cm;
    cm.colors_.clear();
    cm.addColorLast(tgt::Color(0.f, 0.f, 0.f, 1.f));
    cm.addColorLast(tgt::Color(1.f, 0.f, 0.f, 1.f));
    cm.addColorLast(tgt::Color(1.f, 0.8f, 0.f, 1.f));
    cm.setName("Germany");
    return cm;
}

ColorMap ColorMap::createPET() {
    ColorMap cm;
    cm.colors_.clear();
    cm.addColorLast(tgt::Color(0.f, 0.f, 1.f, 1.f));
    cm.addColorLast(tgt::Color(1.f, 0.f, 0.f, 1.f));
    cm.addColorLast(tgt::Color(1.f, 1.f, 0.f, 1.f));
    cm.setName("PET");
    return cm;
}

ColorMap ColorMap::createSequential(const tgt::Color& base) {
    ColorMap cm;
    cm.colors_.clear();
    cm.addColorLast(tgt::Color(1.f, 1.f, 1.f, 1.f));
    cm.addColorLast(base);
    cm.setName("Sequential");
    return cm;
}

ColorMap ColorMap::createFromVector(const std::vector< tgt::Color >& colors) {
    ColorMap cm;
    cm.colors_ = colors;
    return cm;
}

std::vector<std::string> ColorMap::getColorMapLabels() {
    std::vector<std::string> strings;
    strings.push_back(ColorMap::createVisifire1().toString());
    strings.push_back(ColorMap::createVisifire2().toString());
    strings.push_back(ColorMap::createSpectral().toString());
    strings.push_back(ColorMap::createColdHot().toString());
    strings.push_back(ColorMap::createTango().toString());
    strings.push_back(ColorMap::createCandleLight().toString());
    strings.push_back(ColorMap::createSandyShades().toString());
    strings.push_back(ColorMap::createGermany().toString());
    strings.push_back(ColorMap::createPET().toString());
    return strings;
}

ColorMap ColorMap::createColorMap(int index) {
    switch(index) {
        case 0:
            return ColorMap::createVisifire1();
        case 1:
            return ColorMap::createVisifire2();
        case 2:
            return ColorMap::createSpectral();
        case 3:
            return ColorMap::createColdHot();
        case 4:
            return ColorMap::createTango();
        case 5:
            return ColorMap::createCandleLight();
        case 6:
            return ColorMap::createSandyShades();
        case 7:
            return ColorMap::createGermany();
        case 8:
            return ColorMap::createPET();
        default:
            return ColorMap::createSpectral();
    }
}

ColorMap::ColorMap() {
    addColorLast(tgt::Color(0.f, 0.f, 0.f, 1.f));
    addColorLast(tgt::Color(1.f, 1.f, 1.f, 1.f));
}

void ColorMap::serialize(Serializer& s) const {
    s.serialize("colors", colors_);
    s.serialize("name", name_);
}

void ColorMap::deserialize(Deserializer& s) {
    s.deserialize("colors", colors_);
    s.deserialize("name", name_);
}

void ColorMap::addColorLast(const tgt::Color& color) {
    colors_.push_back(color);
}

void ColorMap::addColorAtIndex(const tgt::Color& color, int position) {
    if (position < static_cast<int>(colors_.size()))
        colors_.insert(colors_.begin()+position, color);
}

tgt::Color ColorMap::getColorAtIndex(int position) const {
    position = position % colors_.size();
    return colors_.at(position);
}

tgt::Color ColorMap::getColorAtPosition(float position) const {
    if (position <= 0) {
        return colors_.front();
    }
    else if (position >= 1) {
        return colors_.back();
    }
    else {
        float cPosition = position * (colors_.size()-1);
        tgt::Color left = colors_.at(static_cast<int>(floor(cPosition)));
        tgt::Color right = colors_.at(static_cast<int>(ceil(cPosition)));

        float alpha = static_cast<float>(cPosition - floor(cPosition));
        float iAlpha = 1 - alpha;

        return tgt::Color (iAlpha*left.r + alpha*right.r,
                            iAlpha*left.g + alpha*right.g,
                            iAlpha*left.b + alpha*right.b,
                            iAlpha*left.a + alpha*right.a);
    }
}

int ColorMap::getColorCount() const {
    int toReturn = static_cast<int>(colors_.size());
    return toReturn;
}

ColorMap::InterpolationIterator ColorMap::getInterpolationIterator(int count) const {
    return ColorMap::InterpolationIterator(this, 0, count);
}

ColorMap::GeneratingIterator ColorMap::getGeneratingIterator() const {
    return ColorMap::GeneratingIterator(this, 0);
}

bool ColorMap::operator!=(const ColorMap& rhs) const {
    return !operator==(rhs);
}

bool ColorMap::operator==(const ColorMap& rhs) const {
    if (colors_.size() != rhs.colors_.size())
        return false;

    for (size_t i=0; i <colors_.size(); ++i) {
        if (colors_[i] != rhs.colors_[i])
            return false;
    }
    return true;
}

tgt::vec4 ColorMap::rgbToHsv(tgt::Color rgb){
    float max = std::max(rgb.r, std::max(rgb.g, rgb.b));
    float min = std::min(rgb.r, std::min(rgb.g, rgb.b));

    float h = 0;
    if (max == min)
        h = 0;
    else if (max == rgb.r)
        h = 60*((rgb.g - rgb.b)/(max - min));
    else if (max == rgb.g)
        h = 60*(2 + (rgb.b - rgb.r)/(max - min));
    else
        h = 60*(4 + (rgb.r - rgb.g)/(max - min));

    float s = (max == 0) ? 0 : (max - min)/max;
    float v = max;

    return tgt::vec4(h, s, v, rgb.a);
};

tgt::Color ColorMap::hsvToRgb(tgt::vec4 hsv){
    float hi = floor(hsv.x/60);
    float f = (hsv.x/60) - hi;

    float p = hsv.z * (1-hsv.y);
    float q = hsv.z * (1-(hsv.y*f));
    float t = hsv.z * (1-(hsv.y*(1-f)));

    if (hi == 1)
        return tgt::Color(q, hsv.z, p, hsv.a);
    else if (hi == 2)
        return tgt::Color(p, hsv.z, t, hsv.a);
    else if (hi == 3)
        return tgt::Color(p, q, hsv.z, hsv.a);
    else if (hi == 4)
        return tgt::Color(t, p, hsv.z, hsv.a);
    else if (hi == 5)
        return tgt::Color(hsv.z, p, q, hsv.a);
    else
        return tgt::Color(hsv.z, t, p, hsv.a);
}

std::string ColorMap::toString() const{
    return name_;
}
void ColorMap::setName(const std::string& name) {
    name_ = name;
}

}
