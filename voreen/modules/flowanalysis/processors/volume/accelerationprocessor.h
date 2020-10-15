#ifndef VRN_ACCELERATIONPROCESSOR_H
#define VRN_ACCELERATIONPROCESSOR_H

//add base class header
#include "voreen/core/processors/processor.h"

//add used port headers
#include "voreen/core/ports/volumeport.h"

//add used datastructure headers
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen
{

class VRN_CORE_API AccelerationProcessor : public Processor
{
public:
	AccelerationProcessor();

	virtual Processor *create() const;

	virtual std::string getClassName() const;

	virtual std::string getCategory() const;

	static void Process( const VolumeRAM_Mat3Float& jacobi, const VolumeRAM_3xDouble& velocity, VolumeRAM_3xFloat& outAcceleration );
	static void Process( const VolumeRAM_Mat3Float& jacobi, const VolumeRAM_3xDouble& velocity, VolumeRAM_3xDouble& outAcceleration );

protected:
	virtual void setDescriptions();
	virtual void process();

private:
	VolumePort inportJacobianVolume_;
	VolumePort inportVelocityVolume_;
	VolumePort outport_;
};

} // namespace voreen

#endif // VRN_ACCELERATIONPROCESSOR_H
