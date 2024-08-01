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

#include "wallshearstress.h"

#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "voreen/core/utils/statistics.h"

namespace voreen {

PortDataPointer<GlMeshGeometryUInt32Normal> basicGeometry(const Geometry* geometry) {
    if(auto geom = dynamic_cast<const GlMeshGeometryUInt32Normal*>(geometry)) {
        return PortDataPointer<GlMeshGeometryUInt32Normal>(geom, false);
    }
    else if(auto geom = dynamic_cast<const GlMeshGeometryUInt32NormalTexCoord*>(geometry)) {
        auto g = new GlMeshGeometryUInt32Normal;
        g->setTransformationMatrix(geometry->getTransformationMatrix());
        g->setIndices(geom->getIndices());
        for(const VertexNormalTexCoord& vertex : geom->getVertices()) {
            g->addVertex(vertex.pos_, vertex.normal_);
        }
        return PortDataPointer<GlMeshGeometryUInt32Normal>(g, true);
    }
    return PortDataPointer<GlMeshGeometryUInt32Normal>(nullptr, false);
}

WallShearStress::WallShearStress()
    : Processor()
    , flowInport_(Port::INPORT, "flow", "Flow Input")
    , surfaceInport_(Port::INPORT, "geometry.input", "Geometry Inport")
    , wssSurfaceOutport_(Port::OUTPORT, "geometry.output", "Geometry Outport")
    , intensityVolume_(Port::OUTPORT, "volume.intensities", "Intensity Volume")
    //, dynamicViscosity_("dynamicViscosityProp", "Dynamic viscosity", 0.0035f, 0.f, 0.1f)
    , kinematicViscosity_("kinematicViscosityProp", "Kinematic viscosity (x10^(-6) m^2/s)", 1.35f, 0.1f, 100.0f)
    , density_("densityProp", "Density (kg/m^3)", 988.21f, 0.1f, 10000.0f)
    , maxWSS_("maxWSSProp", "Expected Maximum WSS", 0.00002f, 0.0f, 1000.0f)
    , transferFunction_("transFuncProp", "Color Map")
    , updateVertexPosition_("updateVertexPosition", "Update Vertex Position", true)
    , epsilon_("epsilonProp", "Epsilon (mm)", 0.0f, 0.0f, 100.0f)
    , r_("rProp", "r (mm)", 0.1f, 0.0f, 100.0f)
    , invertNormals_("invertNormalsProp", "Invert normals", false)
    , velocityUnitConversion_("velocityUnitConversionProp", "Input Velocity Unit")
    , interpolation_("interpolation", "Interpolation")
{
    addPort(flowInport_);
    flowInport_.addCondition(new PortConditionVolumeChannelCount(3));
    ON_CHANGE_LAMBDA(flowInport_, [this] {
        if(auto volume = flowInport_.getData()) {
            r_.set(tgt::max(volume->getSpacing()));
        }
    });
    addPort(surfaceInport_);
    addPort(wssSurfaceOutport_);
    addPort(intensityVolume_);

    addProperty(kinematicViscosity_);
    kinematicViscosity_.setNumDecimals(2);

    addProperty(density_);
    density_.setNumDecimals(3);

    addProperty(maxWSS_);
    maxWSS_.setNumDecimals(5);
    addProperty(transferFunction_);

    addProperty(updateVertexPosition_);

    addProperty(epsilon_);
    epsilon_.setNumDecimals(5);
    addProperty(r_);
    r_.setNumDecimals(5);
    addProperty(invertNormals_);

    addProperty(velocityUnitConversion_);
        // Chose the values such that multiplying with real world values we get mm(/s)
        // which (for some reason) is the default voreen length unit.
        velocityUnitConversion_.addOption("km/s", "km/s", 1e6f);
        velocityUnitConversion_.addOption("m/s", "m/s", 1000.0f);
        //velocityUnitConversion_.addOption("dm/s", "dm/s", 100.0f); // Quite unusual.
        velocityUnitConversion_.addOption("cm/s", "cm/s", 10.0f);
        velocityUnitConversion_.addOption("mm/s", "mm/s", 1.0f);
        velocityUnitConversion_.set("m/s");

    addProperty(interpolation_);
    interpolation_.addOption("nearest", "Nearest");
    interpolation_.addOption("linear", "Linear");
    interpolation_.select("linear");
}

Processor* WallShearStress::create() const {
    return new WallShearStress();
}

void WallShearStress::process() {

    // Gather data.
    const VolumeBase* flow = flowInport_.getData();
    VolumeRAMRepresentationLock volumeData(flow);
    tgt::mat4 worldToVoxel = flow->getWorldToVoxelMatrix();
    RealWorldMapping rwm = flow->getRealWorldMapping();

    auto surfaceInput = basicGeometry(surfaceInport_.getData());
    if (!surfaceInput) {
        LERRORC("voreen.WallShearStress", "Surface geometry type not supported");
        return;
    }

    // The intensity volume stores all wss values in a volume to be processed further.
    VolumeRAM_Float* intensityVolume = new VolumeRAM_Float(volumeData->getDimensions());
    intensityVolume->clear();

    const float dynamicViscosity = kinematicViscosity_.get() * 10e-6f * density_.get();
    const float epsilon = epsilon_.get();
    const float r = r_.get();
    const float toMilliMeterPerSecond = velocityUnitConversion_.getValue();
    const bool invertNormals = invertNormals_.get();

    std::function<float(const tgt::vec3&, size_t)> voxel;
    if(interpolation_.get() == "nearest") {
        voxel = [&] (const tgt::vec3& pos, size_t channel) { return volumeData->getVoxelNormalized(pos, channel); };
    }
    else {
        voxel = [&] (const tgt::vec3& pos, size_t channel) { return volumeData->getVoxelNormalizedLinear(pos, channel); };
    }

    auto sample = [&] (tgt::vec3 pos, tgt::vec3 n, float mm = 0.0f) {

        // Transform pos to voxel space.
        pos = worldToVoxel * (pos + n * mm);

        tgt::vec3 v = tgt::vec3::zero;
        for (size_t channel = 0; channel < tgt::vec3::size; channel++) {
            v[channel] = voxel(pos, channel);
            v[channel] = rwm.normalizedToRealWorld(v[channel]);
        }

        return v;
    };

    auto parallel = [] (const tgt::vec3& v, const tgt::vec3& n) {
        return v - (n * tgt::dot(v, n));
    };

    const auto interpolationFromIntensity = [&] ( float intensity ) {

        const auto inputKeys = transferFunction_.get()->getKeys();

        const float div = 1.0f / std::numeric_limits<uint8_t>::max();

        if(intensity <= 0.0f) {
            return tgt::vec4(inputKeys.front()->getColorL()) * div;
        }

        if(intensity >= 1.0f) {
            return tgt::vec4(inputKeys.back()->getColorR()) * div;
        }

        for( size_t i = 0; i < inputKeys.size(); i++ ) {
            if( inputKeys[i]->getIntensity() >= intensity ) {
                float lowerIntensity = 0.0f, upperIntensity = inputKeys[i]->getIntensity();
                tgt::vec4 lowerColor = tgt::vec4::zero, upperColor = inputKeys[i]->getColorL();

                if( i == 0 ) lowerIntensity = 0.0f, lowerColor = inputKeys[0]->getColorR();
                else lowerIntensity = inputKeys[i - 1]->getIntensity(), lowerColor = inputKeys[i - 1]->getColorR();

                float x = ( intensity - lowerIntensity ) / ( upperIntensity - lowerIntensity );

                return (x * upperColor + (1.0f - x) * lowerColor) * div;
            }
        }

        return tgt::vec4::zero;
    };

    const std::vector<VertexNormal>& vertices = surfaceInput->getVertices();
    std::vector<VertexColorNormal> outputVertices(vertices.size());

    std::map<size_t, size_t> access;
    Statistics stats(true);

    // Go through all vertices and calculate magnitude of wss
//#ifdef VRN_MODULE_OPENMP
//#pragma omp parallel for
//#endif
    for (long i=0; i<static_cast<long>(outputVertices.size()); i++) {

        const VertexNormal& vertex = vertices[i];

        tgt::vec3 pos = surfaceInput->getTransformationMatrix() * vertex.pos_;
        tgt::vec3 normal = invertNormals ? tgt::normalize(vertex.normal_) : -tgt::normalize(vertex.normal_);

        tgt::vec3 vectorAtSurface = sample(pos, normal, epsilon);
        tgt::vec3 vectorInLumen = sample(pos, normal, epsilon + r);

        // Calculate WSS magnitude.
        tgt::vec3 v_parallel = parallel(vectorInLumen - vectorAtSurface, normal) * toMilliMeterPerSecond / r;

        float wsr = tgt::length(v_parallel);
        float wss = dynamicViscosity * wsr;

//#pragma omp critical
        {
            tgt::svec3 p = worldToVoxel * pos;

            if(!tgt::hor(tgt::greaterThanEqual(p, intensityVolume->getDimensions()))) {
                size_t idx = VolumeRAM_Float::calcPos(intensityVolume->getDimensions(), p);
                if(intensityVolume->voxel(p) > 0.0f) {

                    access[idx]++;

                    float delta = wss - intensityVolume->voxel(p);
                    intensityVolume->voxel(p) = intensityVolume->voxel(p) + delta / (float)(access[idx] + 1);
                }
                else {
                    intensityVolume->voxel(p) = wss;
                }
            }

            stats.addSample(wss);
        }

        // Normalize to [0, 1].
        wss = (wss / (maxWSS_.get()));

        tgt::vec4 col = interpolationFromIntensity(wss);
        //tgt::vec4 col = tgt::vec4(normal, 1.0f);

        if(updateVertexPosition_.get()) {
            pos = (pos + (normal * epsilon));
        }
        outputVertices[i] = VertexColorNormal(pos, col, vertex.normal_);
    }

    LINFO("Max WSS was " << stats.getMax() << " Pa = " << stats.getMax() * 100.0f << " cPa");
    LINFO("Avg WSS was " << stats.getMean() << " Pa = " << stats.getMean() * 100.0f << " cPa");
    LINFO("Med WSS was " << stats.getMedian() << " Pa = " << stats.getMedian() * 100.0f << " cPa");
    LINFO("Min WSS was " << stats.getMin() << " Pa = " << stats.getMin() * 100.0f << " cPa");

    // Set output.
    GlMeshGeometryUInt32ColorNormal* surfaceOutput = new GlMeshGeometryUInt32ColorNormal();
    surfaceOutput->setVertices(outputVertices);
    surfaceOutput->setIndices(surfaceInput->getIndices());
    wssSurfaceOutport_.setData(surfaceOutput);

    Volume* volume = new Volume(intensityVolume, flow);
    volume->setRealWorldMapping(RealWorldMapping(1.0f, 0.0f, "Pa"));
    intensityVolume_.setData(volume);
}

bool WallShearStress::isReady() const {
    if(!isInitialized()) {
        setNotReadyErrorMessage("Not initialized");
        return false;
    }

    if(!flowInport_.isReady()) {
        setNotReadyErrorMessage("No Flow volume");
        return false;
    }

    if(!surfaceInport_.isReady()) {
        setNotReadyErrorMessage("No surface geometry");
        return false;
    }

    if(!wssSurfaceOutport_.isReady() && !intensityVolume_.isReady()) {
        setNotReadyErrorMessage("No output connected");
        return false;
    }

    return true;
}

}