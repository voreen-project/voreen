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

#include "openlb_parameters.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>

namespace {
    const float PI = 3.14159265358979323846f;
}

namespace voreen {

VelocityCurve::VelocityCurve()
    : periodic_(false)
    , scale_(1.0f)
{
    peakVelocities_[0.0f] = 0.0f;
}

float VelocityCurve::operator()(float t) const {

    if(peakVelocities_.empty()) {
        return 0.0f;
    }

    if(peakVelocities_.size() == 1) {
        return peakVelocities_.begin()->second;
    }

    if (isPeriodic()) {
        t = std::fmod(t - getStartTime(), getEndTime() - getStartTime());
    } else {
        if (t <= peakVelocities_.begin()->first) {
            return peakVelocities_.begin()->second * scale_;
        }

        if (t >= peakVelocities_.rbegin()->first) {
            return peakVelocities_.rbegin()->second * scale_;
        }
    }

    struct Comparator {
        bool operator()(const std::pair<float, float>& p, float value) {
            return p.first < value;
        }
    };

    auto upper = std::lower_bound(peakVelocities_.begin(), peakVelocities_.end(), t, Comparator());
    if (upper == peakVelocities_.begin()) {
        return upper->second * scale_;
    }
    auto lower = upper--;

    const float alpha = (t - lower->first) / (upper->first - lower->first);

    const float value = (1.0f - alpha) * lower->second + alpha * upper->second;

    return value * scale_;
}

float& VelocityCurve::operator[](float t) {
    return peakVelocities_[t];
}

void VelocityCurve::setPeriodic(bool enabled) {
    periodic_ = enabled;
}

bool VelocityCurve::isPeriodic() const {
    return periodic_;
}

void VelocityCurve::setScale(float scale) {
    scale_ = scale;
}

float VelocityCurve::getScale() const {
    return scale_;
}

float VelocityCurve::getMinVelocity() const {
    float min = peakVelocities_.begin()->second;
    for (auto iter = ++peakVelocities_.begin(); iter != peakVelocities_.end(); iter++) {
        min = std::min(iter->second, min);
    }
    return min * scale_;
}

float VelocityCurve::getMaxVelocity() const {
    float max = peakVelocities_.begin()->second;
    for (auto iter = ++peakVelocities_.begin(); iter != peakVelocities_.end(); iter++) {
        max = std::max(iter->second, max);
    }
    return max * scale_;
}

float VelocityCurve::getStartTime() const {
    return peakVelocities_.begin()->first;
}

float VelocityCurve::getEndTime() const {
    return peakVelocities_.rbegin()->first;
}

float VelocityCurve::getDuration() const {
    return peakVelocities_.rbegin()->first - peakVelocities_.begin()->first;
}

VelocityCurve VelocityCurve::createConstantCurve(float value) {
    VelocityCurve curve;
    curve[0.0f] = value;
    return curve;
}

VelocityCurve VelocityCurve::createLinearCurve(float duration, float maxValue) {
    VelocityCurve curve;
    curve[duration] = maxValue;
    return curve;
}

VelocityCurve VelocityCurve::createSinusoidalCurve(float duration, float maxValue, int steps) {
    VelocityCurve curve;
    for(int i=0; i<=steps; i++) {
        float ts = i * duration / steps;
        float value = maxValue * (std::sin(-PI * 0.5f + i * PI / steps) + 1.0f) * 0.5f;
        curve[ts] = value;
    }

    return curve;
}

VelocityCurve VelocityCurve::createHumanHeartBeat() {
    VelocityCurve curve;

    curve[0.040f] =  1.2f / 20.0f * 1.2f;
    curve[0.075f] =  7.0f / 20.0f * 1.2f;
    curve[0.100f] = 17.0f / 20.0f * 1.2f;
    curve[0.125f] = 20.0f / 20.0f * 1.2f;
    curve[0.150f] = 20.0f / 20.0f * 1.2f;
    curve[0.250f] = 10.0f / 20.0f * 1.2f;
    curve[0.280f] =  7.5f / 20.0f * 1.2f;
    curve[0.320f] =  4.2f / 20.0f * 1.2f;
    curve[0.350f] =  1.0f / 20.0f * 1.2f;
    curve[0.380f] =  0.0f / 20.0f * 1.2f;
    curve[0.390f] =  0.0f / 20.0f * 1.2f;
    curve[0.420f] =  0.8f / 20.0f * 1.2f;
    curve[0.455f] =  1.2f / 20.0f * 1.2f;
    curve[0.495f] =  2.5f / 20.0f * 1.2f;
    curve[0.525f] =  2.0f / 20.0f * 1.2f;
    curve[0.560f] =  1.2f / 20.0f * 1.2f;
    curve[0.595f] =  0.3f / 20.0f * 1.2f;
    curve[0.660f] =  0.6f / 20.0f * 1.2f;
    curve[0.700f] =  0.3f / 20.0f * 1.2f;

    return curve;
}

VelocityCurve VelocityCurve::createFromCSV(const std::string& file) {
    VelocityCurve curve;
    curve.peakVelocities_.clear();

    std::ifstream lineStream(file.c_str());
    if (lineStream.fail()) {
        throw std::runtime_error("CSV file could not be opened");
    }

    std::string line;
    std::getline(lineStream, line); // Get (and ignore) header line.
    //std::cout << "Header was: " << line << std::endl;
    while(std::getline(lineStream, line)) {

        std::vector<float> cellValues;

        std::stringstream cellStream(line);
        std::string cell;
        while(std::getline(cellStream, cell, ',')) {
            std::stringstream helper(cell);

            float value = 0.0f;
            helper >> value;
            cellValues.push_back(value);
        }

        if(cellValues.size() != 2) {
            throw std::runtime_error("Invalid tupel size");
        }

        curve[cellValues[0]] = cellValues[1];
    }

    if(curve.peakVelocities_.empty()) {
        throw std::runtime_error("Empty curve");
    }

    return curve;
}

VelocityCurve VelocityCurve::createFromMap(const std::map<float, float>& map) {
    VelocityCurve curve;
    curve.peakVelocities_ = map;
    return curve;
}


Parameters::Parameters()
    : name_()
    , spatialResolution_(0)
    , relaxationTime_(0.0f)
    , characteristicLength_(0.0f)
    , characteristicVelocity_(0.0f)
    , viscosity_(0.0f)
    , density_(0.0f)
    , turbulenceModel_(FTM_NONE)
    , smagorinskyConstant_(0.0f)
    , wallBoundaryCondition_(FBC_NONE)
    , latticePerturbation_(false)
    , inletVelocityMultiplier_(0.0f)
    , geometryIsMesh_(true)
{
}

float Parameters::getReynoldsNumber() const {
    return characteristicVelocity_ * characteristicLength_ / viscosity_;
}

bool Parameters::isValid() const {
    float dx = (characteristicLength_ / spatialResolution_);
    float convertVelocity = 3.0f / (relaxationTime_ - 0.5f) * viscosity_ / dx;
    float uLatticeMax = characteristicVelocity_ * inletVelocityMultiplier_ / convertVelocity;
    return uLatticeMax < 0.4f; // Intrinsic property of LBM.
}

}