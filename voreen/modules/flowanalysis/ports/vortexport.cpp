#include "vortexport.h"
#include <sstream>

namespace voreen
{
	// ------------------- //
	// --- Vortex Port --- //
	// ------------------- //

	VortexPort::VortexPort( PortDirection direction, const std::string& id, const std::string& guiName, bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel )
		: GenericPort<Vortex>( direction, id, guiName, allowMultipleConnections, invalidationLevel )
	{}

	Port* VortexPort::create( PortDirection direction, const std::string& id, const std::string& guiName ) const
	{
		return new VortexPort( direction, id, guiName );
	}
	std::string VortexPort::getClassName() const
	{
		return "VortexPort";
	}
	std::string VortexPort::getContentDescriptionHTML() const
	{
		auto stream = std::stringstream();
		stream << Port::getContentDescriptionHTML();

		if( this->hasData() )
		{
			const auto vortex = this->getData();
			stream << "<br/>Orientation: " << to_string( vortex->getOrientation() );
		}

		return stream.str();
	}

	// ------------------------------ //
	// --- Vortex Collection Port --- //
	// ------------------------------ //

	VortexCollectionPort::VortexCollectionPort( PortDirection direction, const std::string& id, const std::string& guiName, bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel )
		: GenericPort<VortexCollection>( direction, id, guiName, allowMultipleConnections, invalidationLevel )
	{}

	Port* VortexCollectionPort::create( PortDirection direction, const std::string& id, const std::string& guiName ) const
	{
		return new VortexPort( direction, id, guiName );
	}
	std::string VortexCollectionPort::getClassName() const
	{
		return "VortexCollectionPort";
	}
	std::string VortexCollectionPort::getContentDescriptionHTML() const
	{
		auto stream = std::stringstream();
		stream << Port::getContentDescriptionHTML();

		if( this->hasData() )
		{
			const auto collection = this->getData();
			stream << "<br/>Runs: " << std::to_string( collection->runs() );
			stream << "<br/>Timesteps: " << std::to_string( collection->timesteps() );
			stream << "<br/>Vortices: " << std::to_string( collection->totalNumVortices() );
		}

		return stream.str();
	}

	// ------------------------ //
	// --- Vortex List Port --- //
	// ------------------------ //

	VortexListPort::VortexListPort( PortDirection direction, const std::string& id, const std::string& guiName, bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel )
		: GenericPort<std::vector<Vortex>>( direction, id, guiName, allowMultipleConnections, invalidationLevel )
	{}

	Port* VortexListPort::create( PortDirection direction, const std::string& id, const std::string& guiName ) const
	{
		return new VortexListPort( direction, id, guiName );
	}
	std::string VortexListPort::getClassName() const
	{
		return "VortexListPort";
	}
	std::string VortexListPort::getContentDescriptionHTML() const
	{
		auto stream = std::stringstream();
		stream << Port::getContentDescriptionHTML();

		if( this->hasData() )
		{
			const auto vortices = this->getData();
			stream << "<br/>Number of vortices: " << vortices->size();
		}

		return stream.str();
	}
}
