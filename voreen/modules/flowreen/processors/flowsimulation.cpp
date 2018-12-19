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

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"

#include "modules/hdf5/io/hdf5volumereader.h"
#include "modules/hdf5/io/hdf5volumewriter.h"

//NOTE: has to be included at the end due to definitions made in volumeatomic.h!
#include "flowsimulation.h"

#ifndef OLB_PRECOMPILED
#include <olb3D.hh>
#endif

namespace voreen {

/**
 * Helper function to convert a Voreen geometry into an olb STLreader.
 * This is achieved by creating a temporary STL file.
 */
std::unique_ptr<STLreader<T>> convertGeometryToSTL(const Geometry* geometry) {

    if(const GlMeshGeometryBase* data = dynamic_cast<const GlMeshGeometryBase*>(geometry)) {

        if(data->getPrimitiveType() != GL_TRIANGLES) {
            std::cout << "Currently only triangular meshes allowed" << std::endl;
            return nullptr;
        }

        if(!(data->getVertexLayout() & VertexBase::VertexLayout::NORMAL)) {
            std::cout << "Geometry needs to have normals" << std::endl;
            return nullptr;
        }

        // TODO: add this converter to GlMeshGeometry class?
        std::vector<VertexBase> vertices(data->getNumVertices());
        std::vector<uint32_t> indices;
        if(const GlMeshGeometryUInt32Simple* geom = dynamic_cast<const GlMeshGeometryUInt32Simple*>(geometry)) {
            vertices = geom->getVertices();
            indices = geom->getIndices();
        } else if(const GlMeshGeometryUInt32Normal* geom = dynamic_cast<const GlMeshGeometryUInt32Normal*>(geometry)) {
            std::copy(geom->getVertices().begin(), geom->getVertices().end(), vertices.begin());
            indices = geom->getIndices();
        } else if(const GlMeshGeometryUInt32NormalTexCoord* geom = dynamic_cast<const GlMeshGeometryUInt32NormalTexCoord*>(geometry)) {
            std::copy(geom->getVertices().begin(), geom->getVertices().end(), vertices.begin());
            indices = geom->getIndices();
        } else if(const GlMeshGeometryUInt32ColorNormal* geom = dynamic_cast<const GlMeshGeometryUInt32ColorNormal*>(geometry)) {
            std::copy(geom->getVertices().begin(), geom->getVertices().end(), vertices.begin());
            indices = geom->getIndices();
        } else if(const GlMeshGeometryUInt32TexCoord* geom = dynamic_cast<const GlMeshGeometryUInt32TexCoord*>(geometry)) {
            std::copy(geom->getVertices().begin(), geom->getVertices().end(), vertices.begin());
            indices = geom->getIndices();
        } else if(const GlMeshGeometryUInt32ColorNormalTexCoord* geom = dynamic_cast<const GlMeshGeometryUInt32ColorNormalTexCoord*>(geometry)) {
            std::copy(geom->getVertices().begin(), geom->getVertices().end(), vertices.begin());
            indices = geom->getIndices();
        } else {
            std::cout << "Unsupported geometry" << std::endl;
            return nullptr;
        }

        std::string fullName = VoreenApplication::app()->getUniqueTmpFilePath(".stl");
        std::ofstream f(fullName.c_str());
        f << "solid ascii " << fullName << "\n";

        if(data->usesIndexedDrawing()) {
            tgtAssert(data->getNumIndices() % 3 == 0, "No triangle mesh");
            for (size_t i = 0; i < data->getNumIndices() / 3; i+=3) {

                const VertexBase& v0 = vertices[indices[i+0]];
                const VertexBase& v1 = vertices[indices[i+1]];
                const VertexBase& v2 = vertices[indices[i+2]];

                tgt::vec3 normal = tgt::normalize(tgt::cross(v0.pos_-v1.pos_, v0.pos_-v2.pos_));

                f << "facet normal " << normal[0] << " "
                  << normal[1] << " " << normal[2] << "\n";
                f << "    outer loop\n";
                f << "        vertex " << vertices[indices[i+0]].pos_[0] << " "
                  << vertices[indices[i+0]].pos_[1] << " " << vertices[indices[i]].pos_[2]
                  << "\n";
                f << "        vertex " << vertices[indices[i+1]].pos_[0] << " "
                  << vertices[indices[i+1]].pos_[1] << " " << vertices[indices[i]].pos_[2]
                  << "\n";
                f << "        vertex " << vertices[indices[i+2]].pos_[0] << " "
                  << vertices[indices[i+1]].pos_[1] << " " << vertices[indices[i]].pos_[2]
                  << "\n";
                f << "    endloop\n";
                f << "endfacet\n";
            }
        }
        else {
            tgtAssert(data->getNumVertices() % 3 == 0, "No triangle mesh");
            for (size_t i = 0; i < data->getNumVertices() / 3; i+=3) {

                const VertexBase& v0 = vertices[i+0];
                const VertexBase& v1 = vertices[i+1];
                const VertexBase& v2 = vertices[i+2];

                tgt::vec3 normal = tgt::normalize(tgt::cross(v0.pos_-v1.pos_, v0.pos_-v2.pos_));

                f << "facet normal " << normal[0] << " "
                  << normal[1] << " " << normal[2] << "\n";
                f << "    outer loop\n";
                f << "        vertex " << vertices[i+0].pos_[0] << " "
                  << vertices[i+0].pos_[1] << " " << vertices[i].pos_[2]
                  << "\n";
                f << "        vertex " << vertices[i+1].pos_[0] << " "
                  << vertices[i+1].pos_[1] << " " << vertices[i].pos_[2]
                  << "\n";
                f << "        vertex " << vertices[i+2].pos_[0] << " "
                  << vertices[i+1].pos_[1] << " " << vertices[i].pos_[2]
                  << "\n";
                f << "    endloop\n";
                f << "endfacet\n";
            }
        }

        f.close();

        std::unique_ptr<STLreader<T>> reader;
        try {
            reader.reset(new STLreader<T>(fullName, 1.0f));
        }
        catch(const std::runtime_error& error) {
            std::cout << error.what() << std::endl;
        }
        return reader;
    }
    else {
        std::cout << "Geometry not supported!" << std::endl;
        return nullptr;
    }
}


const std::string FlowSimulation::loggerCat_("voreen.flowreen.FlowSimulation");

FlowSimulation::FlowSimulation()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    // ports
    , geometryDataPort_(Port::INPORT, "geometryDataPort", "Geometry Input", false)
    , measuredDataPort_(Port::INPORT, "measuredDataPort", "Measured Data Input", false)
    , outport_(Port::OUTPORT, "outport", "Time Series Output")
    , simulationTime_("simulationTime", "Simulation Time (s)", 2.0f, 0.1f, 10.0f)
    , temporalResolution_("temporalResolution", "Temporal Resolution (ms)", 3.1f, 1.0f, 30.0f)
    , characteristicLength_("characteristicLength", "Characteristic Length (mm)", 22.46f, 1.0f, 100.0f)
    , viscosity_("viscosity", "Viscosity (m^2/s)", 3.5e-10, 3e-10, 4e-10)
    , density_("density", "Density (kg/m^3)", 1000.0f, 1000.0f, 1100.0f)
{
    addPort(geometryDataPort_);
    addPort(measuredDataPort_);
    addPort(outport_);

    addProperty(simulationTime_);
    addProperty(temporalResolution_);
    addProperty(characteristicLength_);
    addProperty(viscosity_);
    addProperty(density_);
}

FlowSimulation::~FlowSimulation() {
}

bool FlowSimulation::isReady() const {
    if(!isInitialized()) {
        setNotReadyErrorMessage("Not initialized.");
        return false;
    }
    if(!geometryDataPort_.isReady()) {
        setNotReadyErrorMessage("Geometry Port not ready.");
        return false;
    }
    // Note: measuredDataPort ist optional!
    return true;
}

FlowSimulationInput FlowSimulation::prepareComputeInput() {
    const Geometry* geometryData = geometryDataPort_.getData();
    if (!geometryData) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    const VolumeList* measuredData = measuredDataPort_.getData();
    if(!measuredData || measuredData->empty()) {
        throw InvalidInputException("Unsteered simulations currently not supported", InvalidInputException::S_ERROR);
    }

    LINFO("Configuring a steered simulation");

    // TODO: create new data-/port- type which resamples all contained volume into a cube or at least performs the checks below.
    // Check for volume compatibility
    VolumeBase* volumeT0 = measuredData->first();
    // Currently only 3xFloat Volumes are considered. This condition could be relaxed in the future.
    if(volumeT0->getFormat() != VolumeGenerator3xFloat().getFormat()) {
        throw InvalidInputException("Measured data contains volume different from 3xFloat", InvalidInputException::S_ERROR);
    }
    if(volumeT0->getDimensions() != tgt::svec3(volumeT0->getDimensions().x)) {
        throw InvalidInputException("Measured data must have dimensions: n x n x n", InvalidInputException::S_ERROR);
    }
    if(volumeT0->getSpacing() != tgt::vec3(volumeT0->getSpacing().x)) {
        throw InvalidInputException("Measured data must have spacing: n x n x n", InvalidInputException::S_ERROR);
    }

    if(!volumeT0->hasDerivedData<VolumeMinMaxMagnitude>()) {
        LWARNING("Calculating VolumeMinMaxMagnitude. This may take a while...");
        //throw InvalidInputException("VolumeMinMaxMagnitude not available!", InvalidInputException::S_WARNING);
    }

    float minVelocityMagnitude = volumeT0->getDerivedData<VolumeMinMaxMagnitude>()->getMinMagnitude();
    float maxVelocityMagnitude = volumeT0->getDerivedData<VolumeMinMaxMagnitude>()->getMaxMagnitude();

    for(size_t i=1; i<measuredData->size(); i++) {
        VolumeBase* volumeTi = measuredData->at(i);
        if(volumeT0->getFormat() != volumeTi->getFormat()
            || volumeT0->getDimensions() != volumeTi->getDimensions()
            || volumeT0->getSpacing() != volumeTi->getSpacing()) {
            throw InvalidInputException("Measured data contains different kinds of volumes.", InvalidInputException::S_ERROR);
        }

        minVelocityMagnitude = std::min(minVelocityMagnitude, volumeTi->getDerivedData<VolumeMinMaxMagnitude>()->getMinMagnitude());
        maxVelocityMagnitude = std::min(maxVelocityMagnitude, volumeTi->getDerivedData<VolumeMinMaxMagnitude>()->getMaxMagnitude());
    }

    // === 1st Step: Initialization ===
    UnitConverter<T,DESCRIPTOR> converter(
            (T)   volumeT0->getSpacing().x,             // physDeltaX: spacing between two lattice cells in __m__
            (T)   temporalResolution_.get(),            // physDeltaT: time step in __s__
            (T)   characteristicLength_.get(),          // charPhysLength: reference length of simulation geometry
            (T)   maxVelocityMagnitude/1000.0,          // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
            (T)   viscosity_.get(),                     // physViscosity: physical kinematic viscosity in __m^2 / s__
            (T)   density_.get()                        // physDensity: physical density in __kg / m^3__
    );
    // Prints the converter log as console output
    converter.print();
    // Writes the converter log in a file
    converter.write("aorta3d");

    // === 2nd Step: Prepare Geometry ===

    // Instantiation of the STLreader class
    // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.
    std::unique_ptr<STLreader<T>> stlReader = convertGeometryToSTL(geometryData);
    if(!stlReader) {
        throw InvalidInputException("Geometry could not be initialized", InvalidInputException::S_ERROR);
    }
    //STLreader<T> stlReader( geometryFile_.get(), converter.getConversionFactorLength(), 0.001, 0, true );
/*
    const std::string volumeLocation = HDF5VolumeWriter::VOLUME_DATASET_NAME;
    const std::string baseType = volumeT0->getBaseType();
    const tgt::svec3 outputDim = volumeT0->getDimensions();
    const int deflateLevel = 1;

    std::vector<std::unique_ptr<HDF5FileVolume>> outputVolumes;
    for(size_t i = 0; i<10; i++) {
        std::unique_ptr<HDF5FileVolume> outputVolume = nullptr;
        try {
            outputVolume = std::unique_ptr<HDF5FileVolume>(
                    HDF5FileVolume::createVolume(volumeFilePath, volumeLocation, baseType, outputDim, 1, true,
                                                 deflateLevel, tgt::svec3(outputDim.xy(), 1), false));
        } catch (tgt::IOException e) {
            throw InvalidInputException("Could not create output volume.", InvalidInputException::S_ERROR);
        }

        outputVolume->writeSpacing(volumeT0->getSpacing());
        outputVolume->writeOffset(volumeT0->getOffset());
        outputVolume->writePhysicalToWorldTransformation(volumeT0->getPhysicalToWorldMatrix());
        outputVolume->writeRealWorldMapping(volumeT0->getRealWorldMapping());
        outputVolumes.emplace_back(std::move(outputDim));
    }
*/
    return FlowSimulationInput{
            simulationTime_.get(),
            converter,
            std::move(stlReader),
            //outputVolumes
    };
}

FlowSimulationOutput FlowSimulation::compute(FlowSimulationInput input, ProgressReporter& progressReporter) const {

    progressReporter.setProgress(0.0f);

    std::unique_ptr<VolumeList> output = nullptr;//std::move(input.outputVolumes);

    // Needs to be initialized in each new thread to be used.
    olb::olbInit(nullptr, nullptr);

    IndicatorLayer3D<T> extendedDomain( *input.stlReader, input.converter.getConversionFactorLength() );

    // Instantiation of a cuboidGeometry with weights
    const int noOfCuboids = 2;
    CuboidGeometry3D<T> cuboidGeometry( extendedDomain, input.converter.getConversionFactorLength(), noOfCuboids );

    // Instantiation of a loadBalancer
    HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

    // Instantiation of a superGeometry
    SuperGeometry3D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );

    prepareGeometry( input.converter, extendedDomain, *input.stlReader, superGeometry );

    // === 3rd Step: Prepare Lattice ===
    SuperLattice3D<T, DESCRIPTOR> sLattice(superGeometry);

    SmagorinskyBGKdynamics<T, DESCRIPTOR> bulkDynamics(
            input.converter.getLatticeRelaxationFrequency(),
            instances::getBulkMomenta<T, DESCRIPTOR>(), 0.1 );

    // choose between local and non-local boundary condition
    sOnLatticeBoundaryCondition3D<T,DESCRIPTOR> sBoundaryCondition( sLattice );
    createInterpBoundaryCondition3D<T,DESCRIPTOR>( sBoundaryCondition );
    // createLocalBoundaryCondition3D<T,DESCRIPTOR>(sBoundaryCondition);

    sOffLatticeBoundaryCondition3D<T, DESCRIPTOR> sOffBoundaryCondition( sLattice );
    createBouzidiBoundaryCondition3D<T, DESCRIPTOR> ( sOffBoundaryCondition );

    prepareLattice( sLattice, input.converter, bulkDynamics,
                    sBoundaryCondition, sOffBoundaryCondition,
                    *input.stlReader, superGeometry );

    // === 4th Step: Main Loop ===
    for ( int iT = 0; iT <= input.converter.getLatticeTime( input.simulationTime ); iT++ ) {

        // === 5th Step: Definition of Initial and Boundary Conditions ===
        setBoundaryValues( sLattice, sOffBoundaryCondition, input.converter, iT, superGeometry );

        // === 6th Step: Collide and Stream Execution ===
        sLattice.collideAndStream();

        // === 7th Step: Computation and Output of the Results ===
        bool success = getResults(sLattice, input.converter, iT, output.get());
        if(!success)
            break;

        float progress = iT / (input.converter.getLatticeTime( input.simulationTime ) + 1.0f);
        progressReporter.setProgress(progress);
    }
    progressReporter.setProgress(1.0f);

    // Done.
    return FlowSimulationOutput{
            std::move(output)
    };
}

void FlowSimulation::processComputeOutput(FlowSimulationOutput output) {
}

// Stores data from stl file in geometry in form of material numbers
void FlowSimulation::prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter, IndicatorF3D<T>& indicator,
                                      STLreader<T>& stlReader, SuperGeometry3D<T>& superGeometry ) const {

    LINFO("Prepare Geometry ...");

    // TODO

    LINFO("Prepare Geometry ... OK");
}

// Set up the geometry of the simulation
void FlowSimulation::prepareLattice( SuperLattice3D<T, DESCRIPTOR>& lattice,
                                     UnitConverter<T,DESCRIPTOR> const& converter,
                                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                                     sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>& bc,
                                     sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& offBc,
                                     STLreader<T>& stlReader, SuperGeometry3D<T>& superGeometry ) const {

    LINFO("Prepare Lattice ...");

    // TODO

    LINFO("Prepare Lattice ... OK");
}

// Generates a slowly increasing sinuidal inflow
void FlowSimulation::setBoundaryValues( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                        sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& offBc,
                                        UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                                        SuperGeometry3D<T>& superGeometry ) const {
    // TODO
}

// Computes flux at inflow and outflow
bool FlowSimulation::getResults( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                 UnitConverter<T,DESCRIPTOR>& converter, int iT,
                                 VolumeList* volumeList ) const {
    OstreamManager clout( std::cout,"getResults" );

    SuperVTMwriter3D<T> vtmWriter( "aorta3d" );
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );

    const int vtkIter  = converter.getLatticeTime( .1 );
    const int statIter = converter.getLatticeTime( .1 );

    // Writes the vtk files
    if ( iT%vtkIter==0 ) {
        // TODO
        //sLattice.get(0, 0, 0, 0).computeU();
    }

    /*
    if ( sLattice.getStatistics().getMaxU() > 0.3 ) {
        clout << "PROBLEM uMax=" << sLattice.getStatistics().getMaxU() << std::endl;
        vtmWriter.write( iT );
        std::exit( 0 );
    }
    */
}

}   // namespace
