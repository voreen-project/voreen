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

#include "serializationhelper.h"

namespace voreen {

VelocityCurveSerializable::VelocityCurveSerializable() {
}

VelocityCurveSerializable::VelocityCurveSerializable(const VelocityCurve& other)
    : VelocityCurve(other)
{
}

void VelocityCurveSerializable::serialize(Serializer& s) const {
    s.serialize("peakVelocities", peakVelocities_);
    s.serialize("periodic", periodic_);
    s.serialize("scale", scale_);
}

void VelocityCurveSerializable::deserialize(Deserializer& s) {
    s.deserialize("peakVelocities", peakVelocities_);
    s.deserialize("periodic", periodic_);
    s.optionalDeserialize("scale", scale_, 1.0f); // Maintain backwards compatibility.
}

FlowIndicatorSerializable::FlowIndicatorSerializable() {
}

FlowIndicatorSerializable::FlowIndicatorSerializable(const FlowIndicator& other)
    : FlowIndicator(other)
{
}

void FlowIndicatorSerializable::serialize(Serializer& s) const {
    s.serialize("type_", type_); // "type" is a reserved xml keyword.
    s.serialize("id_", id_); // "id" is a reserved xml keyword.
    s.serialize("name", name_);
    s.serialize("color", color_);
    s.serialize("center", center_);
    s.serialize("normal", normal_);
    s.serialize("radius", radius_);
    s.serialize("length", length_);
    s.serialize("flowProfile", flowProfile_);
    s.serialize("velocityCurve", VelocityCurveSerializable(velocityCurve_));
}

void FlowIndicatorSerializable::deserialize(Deserializer& s) {
    int type = FIT_CANDIDATE;
    s.deserialize("type_", type);
    type_ = static_cast<FlowIndicatorType>(type);
    s.deserialize("id_", id_);
    s.optionalDeserialize("name", name_, "Indicator " + std::to_string(id_));
    s.optionalDeserialize("color", color_, tgt::vec4(1.0));
    s.deserialize("center", center_);
    s.deserialize("normal", normal_);
    s.deserialize("radius", radius_);
    s.optionalDeserialize("length", length_, 0.0f);
    int profile = FP_POISEUILLE;
    s.deserialize("flowProfile", profile);
    flowProfile_ = static_cast<FlowProfile>(profile);
    VelocityCurveSerializable curve;
    s.deserialize("velocityCurve", curve);
    velocityCurve_ = static_cast<VelocityCurve>(curve);
}

ParametersSerializable::ParametersSerializable() {
}

ParametersSerializable::ParametersSerializable(const Parameters& other)
        : Parameters(other) {
}

void ParametersSerializable::serialize(Serializer& s) const {
    s.serialize("name", name_);
    s.serialize("spatialResolution", spatialResolution_);
    s.serialize("relaxationTime", relaxationTime_);
    s.serialize("characteristicLength", characteristicLength_);
    s.serialize("characteristicVelocity", characteristicVelocity_);
    s.serialize("viscosity", viscosity_);
    s.serialize("density", density_);
    s.serialize("turbulenceModel", turbulenceModel_);
    s.serialize("smagorinskyConstant", smagorinskyConstant_);
    s.serialize("wallBoundaryCondition", wallBoundaryCondition_);
    s.serialize("latticePerturbation", latticePerturbation_);
    s.serialize("inletVelocityMultiplier", inletVelocityMultiplier_);
}

void ParametersSerializable::deserialize(Deserializer& s) {
    s.deserialize("name", name_);
    s.deserialize("spatialResolution", spatialResolution_);
    s.deserialize("relaxationTime", relaxationTime_);
    s.deserialize("characteristicLength", characteristicLength_);
    s.deserialize("characteristicVelocity", characteristicVelocity_);
    s.deserialize("viscosity", viscosity_);
    s.deserialize("density", density_);
    int turbulenceModel = FTM_NONE;
    s.optionalDeserialize("turbulenceModel", turbulenceModel, static_cast<int>(FTM_SMAGORINSKY));
    turbulenceModel_ = static_cast<FlowTurbulenceModel>(turbulenceModel);
    s.deserialize("smagorinskyConstant", smagorinskyConstant_);
    int wallBoundaryCondition = FBC_NONE;
    s.optionalDeserialize("wallBoundaryCondition", wallBoundaryCondition, static_cast<int>(FBC_BOUZIDI));
    wallBoundaryCondition_ = static_cast<FlowBoundaryCondition>(wallBoundaryCondition);
    s.optionalDeserialize("latticePerturbation", latticePerturbation_, false);
    s.optionalDeserialize("inletVelocityMultiplier", inletVelocityMultiplier_, 1.0f);
}

}