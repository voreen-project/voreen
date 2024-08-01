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

#include "pyrogenicvorticitymapper.h"

namespace voreen {

const std::string PyrogenicVorticityMapper::loggerCat_ = "PyrogenicVorticityMapper";

struct AlgorithmInput {
    tgt::vec3 vorticity;
    tgt::vec3 materialDerivative;
    bool hasMaterialDerivative;
    float materialDerivativeEpsilon;
    float threshold;
};

PyrogenicVorticityMapper::PyrogenicVorticityMapper() 
    : Processor()
    , vorticityInport_( Port::INPORT, "jacobianInport", "Jacobian of the volume of interest" )
    , materialDerivative_(Port::INPORT, "materialDerivative", "Optional material derivative input")
    , outputVolume_( Port::OUTPORT, "outputVolume", "Scalar representation of the pyrogenic vorticity field." )
    , mappingAlgorithm_("mappingAlgorithm", "Mapping Algorithm: ", Processor::InvalidationLevel::INVALID_PARAMETERS)
    , pyroVorticityVolume_(Port::OUTPORT, "pyrogenicVorticity", "Pyrogenic part of the vorticity field")
    , enableFirstComponent_("enableFirstComponent", "First Component: ", true)
    , enableSecondComponent_("enableSecondComponent", "Second Component: ", false)
    , enableThirdComponent_("enableThirdComponent", "Third Component: ", true)
    , epsilonProperty_("epsilon", "Epsilon", 0.f, 0.f, 10.f)
    , invertProperty_("invertOutout", "Invert Output: ")
{

    addPort(vorticityInport_);
    addPort(materialDerivative_);
    addPort(outputVolume_);
    addPort(pyroVorticityVolume_);

    addProperty(mappingAlgorithm_);

    addProperty(enableFirstComponent_);
    addProperty(enableSecondComponent_);
    addProperty(enableThirdComponent_);

    addProperty(epsilonProperty_);
    addProperty(invertProperty_);

    mappingAlgorithm_.addOption("length", "Map to magnitude", 0);
    mappingAlgorithm_.addOption("pyrogenic * length", "Map to pyrogenic value by magnitude", 1);
    mappingAlgorithm_.addOption("pyrogenic", "Map to pyrogenic value", 2);
    mappingAlgorithm_.addOption("sgn(x)", "Map to the sign of the x component", 3);
    mappingAlgorithm_.addOption("sgn(y)", "Map to the sign of the y component", 4);
    mappingAlgorithm_.addOption("sgn(z)", "Map to the sign of the z component", 5);
    mappingAlgorithm_.addOption("sgn(x) * length", "Map to the sign of the x component scaled by the magnitude", 6);
    mappingAlgorithm_.addOption("sgn(y) * length", "Map to the sign of the y component scaled by the magnitude", 7);
    mappingAlgorithm_.addOption("sgn(z) * length", "Map to the sign of the z component scaled by the magnitude", 8);
}

bool PyrogenicVorticityMapper::isReady() const {
    return 
        vorticityInport_.hasData() &&
        outputVolume_.isConnected();
}

float map_magnitude(const AlgorithmInput& input) {
    return std::abs(tgt::lengthSq(input.vorticity));
}

float map_pyrogenicBinary(const AlgorithmInput& input) {
    // test for
    // |Dw/Dt| > epsilon 
    if (input.hasMaterialDerivative && std::abs(input.materialDerivative.y) > input.materialDerivativeEpsilon)
        return 0.f;
    
    if (input.vorticity.x > input.threshold && input.vorticity.z > input.threshold)
        return 1;
    else if (input.vorticity.x < -input.threshold && input.vorticity.z < -input.threshold)
        return -1;
    
    return 0.f;
}

float map_pyrogenic(const AlgorithmInput& input) {
    auto magnitude = map_magnitude(input);
    return magnitude * map_pyrogenicBinary(input);
}

float map_signXBinary(const AlgorithmInput& input) {
    if (input.vorticity.x < -input.threshold) return -1.f;
    if (input.vorticity.x > input.threshold) return 1.f;
    return 0.f;
}

float map_signYBinary(const AlgorithmInput& input) {
    if (input.vorticity.y < -input.threshold) return -1.f;
    if (input.vorticity.y > input.threshold) return 1.f;
    return 0.f;
}

float map_signZBinary(const AlgorithmInput& input) {
    if (input.vorticity.z < -input.threshold) return -1.f;
    if (input.vorticity.z > input.threshold) return 1.f;
    return 0.f;
}

float map_signX(const AlgorithmInput& input) {
    return map_signXBinary(input) * map_magnitude(input);
}

float map_signY(const AlgorithmInput& input) {
    return map_signYBinary(input) * map_magnitude(input);
}

float map_signZ(const AlgorithmInput& input) {
    return map_signZBinary(input) * map_magnitude(input);
}

void PyrogenicVorticityMapper::process() {
    auto inputVolume = vorticityInport_.getData();
    VolumeRAMRepresentationLock vorticityVolume(inputVolume);
    const auto* vorticityVolumeData = dynamic_cast<const VolumeRAM_3xFloat*>(*vorticityVolume);
    if (vorticityVolumeData == nullptr)
        throw std::runtime_error("Expected vorticity volume as inport!");

    const VolumeRAM_3xFloat* matDerivData = nullptr;

    if (materialDerivative_.hasData()) {
        auto matDerivPortData = materialDerivative_.getData();
        VolumeRAMRepresentationLock matDerivVolume(matDerivPortData);
        matDerivData = dynamic_cast<const VolumeRAM_3xFloat*>(*matDerivVolume);
    }

    tgt::Vector3<long> dimensions = inputVolume->getDimensions();
    RealWorldMapping rwm = inputVolume->getRealWorldMapping();

    auto mappingVolume = new VolumeRAM_Float(inputVolume->getDimensions());
    auto pVorticityVolume = new VolumeRAM_3xFloat(inputVolume->getDimensions());

    auto algorithms = std::vector<float(*)(const AlgorithmInput&)>{
        map_magnitude,
        map_pyrogenic,
        map_pyrogenicBinary,
        map_signXBinary,
        map_signYBinary,
        map_signZBinary,
        map_signX,
        map_signY,
        map_signZ
    };

    AlgorithmInput input;
    input.hasMaterialDerivative = matDerivData != nullptr;
    input.threshold = epsilonProperty_.get();
    input.materialDerivativeEpsilon = epsilonProperty_.get();

    auto algorithm = algorithms[ mappingAlgorithm_.getValue() ];

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
    for (long z = 0; z < dimensions.z; z++) {
        for (long y = 0; y < dimensions.y; y++) {
            for (long x = 0; x < dimensions.x; x++) {
                tgt::svec3 pos(x, y, z);
                auto vorticity = vorticityVolumeData->voxel(pos);

                auto pyrogenicVorticity = tgt::vec3{
                    vorticity.x * enableFirstComponent_.get(), 
                    vorticity.y * enableSecondComponent_.get(),  
                    vorticity.z * enableThirdComponent_.get()
                };

                float materialDerivative = 0.f;
                if (matDerivData)
                    input.materialDerivative = matDerivData->voxel(pos);
                
                input.vorticity = pyrogenicVorticity;

                auto mapped = algorithm(input);

                if (invertProperty_.get()) {
                    mappingVolume->voxel(pos) = -mapped;
                }
                else {
                    mappingVolume->voxel(pos) = mapped;
                }
                pVorticityVolume->voxel(pos) = pyrogenicVorticity;
            }
        }
    }

    auto* volume1 = new Volume(mappingVolume, inputVolume);
    volume1->setRealWorldMapping(RealWorldMapping()); // Override to default rwm.
    volume1->setModality(Modality("pyrogenic_mapping"));
    outputVolume_.setData(volume1);
    if (pyroVorticityVolume_.isConnected()) {
        auto* volume2 = new Volume(pVorticityVolume, inputVolume);
        volume2->setRealWorldMapping(RealWorldMapping()); // Override to default rwm.
        volume2->setModality(Modality("pyrogenic_mapping"));
        pyroVorticityVolume_.setData(volume2);
    }
}

}
