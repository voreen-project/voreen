#pragma once
#include <string>
#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"

#include "modules/flowanalysis/ports/parallelvectorsolutionpointsport.h"

namespace voreen
{
	class ParallelVectors : public Processor
	{
	public:
		ParallelVectors();
		virtual Processor *create() const { return new ParallelVectors(); }
		virtual std::string getClassName() const { return "ParallelVectors"; }
		virtual std::string getCategory() const { return "Volume Processing"; }
		virtual bool isReady() const;

		static void Process( const VolumeRAM_3xFloat& V, const VolumeRAM_3xFloat& W, const VolumeRAM_Mat3Float* jacobi, const VolumeRAM* mask, ParallelVectorSolutions& outSolution );
		static void Process( const VolumeRAM_3xDouble& V, const VolumeRAM_3xDouble& W, const VolumeRAM_Mat3Float* jacobi, const VolumeRAM* mask, ParallelVectorSolutions& outSolution );

		static constexpr auto TetrahedraPerCube = 6;
		static constexpr auto TrianglesPerTetrahedron = 4;

	protected:
		virtual void process();

	private:
		void onChangedJacobianData();
		VolumePort _inV, _inW, _inJacobi, _inMask;
		ParallelVectorSolutionPointsPort _out;
		BoolProperty _sujudiHaimes;

		using Triangle = std::array<tgt::svec3, 3>;
		using Tet = std::array<Triangle, TrianglesPerTetrahedron>;
	};
} // namespace voreen
