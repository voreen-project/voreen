#include "accelerationprocessor.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include <chrono>

namespace voreen
{
AccelerationProcessor::AccelerationProcessor() : Processor(), inportJacobianVolume_(Port::INPORT, "inportJacobianVolume", "3x3 Float Matrix Volume (Jacobian)"), inportVelocityVolume_(Port::INPORT, "inportVelocityVolume", "3-Double Vector Volume (Velocity)"), outport_(Port::OUTPORT, "outport", "3-Float Vector Volume (Acceleration)")
{
	inportJacobianVolume_.addCondition(new PortConditionVolumeType("Matrix3(float)", "Volume_Mat3Float"));
	inportVelocityVolume_.addCondition(new PortConditionVolumeType3xDouble);
	addPort(inportJacobianVolume_);
	addPort(inportVelocityVolume_);
	addPort(outport_);
}

Processor *AccelerationProcessor::create() const
{
	return new AccelerationProcessor();
}

std::string AccelerationProcessor::getClassName() const
{
	return "AccelerationProcessor";
}

std::string AccelerationProcessor::getCategory() const
{
	return "Volume Processing";
}

void AccelerationProcessor::Process( const VolumeRAM_Mat3Float& jacobi, const VolumeRAM_3xDouble& velocity, VolumeRAM_3xFloat& outAcceleration )
{
	const auto timeBegin = std::chrono::high_resolution_clock::now();

#pragma omp parallel for
	for( long i = 0; i < static_cast<long>( jacobi.getNumVoxels() ); ++i )
	{
		outAcceleration.voxel( i ) = jacobi.voxel( i ) * velocity.voxel( i );
	}

	const auto timeEnd = std::chrono::high_resolution_clock::now();
	const auto time = std::chrono::duration_cast<std::chrono::milliseconds>( timeEnd - timeBegin ).count();
	// std::cout << "[AccelerationProcessor]: Finished in " << time << " ms." << std::endl;
}
void AccelerationProcessor::Process( const VolumeRAM_Mat3Float& jacobi, const VolumeRAM_3xDouble& velocity, VolumeRAM_3xDouble& outAcceleration )
{
	const auto timeBegin = std::chrono::high_resolution_clock::now();

#pragma omp parallel for
	for( long i = 0; i < static_cast<long>( jacobi.getNumVoxels() ); ++i )
	{
		outAcceleration.voxel( i ) = jacobi.voxel( i ) * velocity.voxel( i );
	}

	const auto timeEnd = std::chrono::high_resolution_clock::now();
	const auto time = std::chrono::duration_cast<std::chrono::milliseconds>( timeEnd - timeBegin ).count();
	// std::cout << "[AccelerationProcessor]: Finished in " << time << " ms." << std::endl;
}

void AccelerationProcessor::setDescriptions()
{
	setDescription("Computes acceleration volume by multiplying jacobian and velocity for each voxel");
}

void AccelerationProcessor::process()
{
	auto jacobianVolume = dynamic_cast<const VolumeRAM_Mat3Float *>(inportJacobianVolume_.getData()->getRepresentation<VolumeRAM>());
	auto velocityVolume = dynamic_cast<const VolumeRAM_3xDouble *>(inportVelocityVolume_.getData()->getRepresentation<VolumeRAM>());

	auto accelerationVolume = std::unique_ptr<VolumeRAM_3xDouble>( new VolumeRAM_3xDouble(jacobianVolume->getDimensions()) );

	AccelerationProcessor::Process( *jacobianVolume, *velocityVolume, *accelerationVolume );
	
	Volume* output = new Volume(accelerationVolume.release(), inportVelocityVolume_.getData()->getSpacing(), inportVelocityVolume_.getData()->getOffset());
	output->getMetaDataContainer().addMetaData("name", new StringMetaData("Acceleration"));
	outport_.setData(output);

}

} // namespace voreen
