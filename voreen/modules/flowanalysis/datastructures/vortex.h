#pragma once
#ifndef VRN_VORTEX_H
#define VRN_VORTEX_H

#include "tgt/vector.h"
#include <vector>

namespace voreen
{
	class Vortex
	{
	public:
		enum class Orientation : int32_t
		{
			eUnknown = 0x0,
			eClockwise = 0x1,
			eCounterClockwise = 0x2
		};

		Vortex() = default;
		Vortex( Orientation orientation, std::vector<tgt::vec3> coreline );
		Vortex( std::istream& stream );

		void serialize( std::ostream& stream ) const;

		void setOrientation( Orientation orientation ) noexcept;
		Orientation getOrientation() const noexcept;

		void setCoreline( std::vector<tgt::vec3> coreline ) noexcept;
		const std::vector<tgt::vec3>& coreline() const noexcept;

	private:
		Orientation _orientation;
		std::vector<tgt::vec3> _coreline;
	};

	class VortexCollection
	{
	public:
		struct VortexID
		{
			size_t run, timestep, index;

			VortexID() noexcept = default;
			VortexID( size_t run, size_t timestep, size_t index );

			bool operator==( const VortexID& other ) const noexcept;
			bool operator!=( const VortexID& other ) const noexcept;

			static VortexID Invalid;
		};

		VortexCollection() = default;
		VortexCollection( size_t runs, size_t timesteps );
		VortexCollection( std::istream& stream );

		void serialize( std::ostream& stream ) const;

		size_t runs() const noexcept;
		size_t timesteps() const noexcept;
		size_t totalNumVortices() const;

		const std::vector<Vortex>& vortices( size_t run, size_t timestep ) const;
		void setVortices( size_t run, size_t timestep, std::vector<Vortex> vortices );

		const std::vector<VortexID>& matches( size_t run, size_t timestep, size_t index ) const;
		const std::vector<VortexID>& matches( VortexID vortexID ) const;

		void addMatch( VortexID first, VortexID second );

	private:
		size_t _runs, _timesteps;
		std::vector<std::vector<Vortex>> _vortices;
		std::vector<std::vector<std::vector<VortexID>>> _matches;
	};

	std::string to_string( Vortex::Orientation orientation );
}

#endif // VRN_VORTEX_H