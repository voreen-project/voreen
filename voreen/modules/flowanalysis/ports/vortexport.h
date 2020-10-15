#pragma once
#ifndef VRN_VORTEXPORT_H
#define VRN_VORTEXPORT_H

#include "voreen/core/ports/genericport.h"
#include "../datastructures/vortex.h"

namespace voreen
{
#ifdef DLL_TEMPLATE_INST
	template class VRN_CORE_API GenericPort<Vortex>;
#endif

	class VRN_CORE_API VortexPort : public GenericPort<Vortex>
	{
	public:
		VortexPort( PortDirection direction, const std::string& id, const std::string& guiName = {}, bool allowMultipleConnections = false, Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT );

		virtual Port* create( PortDirection direction, const std::string& id, const std::string& guiName = {} ) const;
		virtual std::string getClassName() const;
		virtual std::string getContentDescriptionHTML() const;
	};

#ifdef DLL_TEMPLATE_INST
	template class VRN_CORE_API GenericPort<VortexCollection>;
#endif

	class VRN_CORE_API VortexCollectionPort : public GenericPort<VortexCollection>
	{
	public:
		VortexCollectionPort( PortDirection direction, const std::string& id, const std::string& guiName = {}, bool allowMultipleConnections = false, Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT );

		virtual Port* create( PortDirection direction, const std::string& id, const std::string& guiName = {} ) const;
		virtual std::string getClassName() const;
		virtual std::string getContentDescriptionHTML() const;
	};

#ifdef DLL_TEMPLATE_INST
	template class VRN_CORE_API GenericPort<std::vector<Vortex>>;
#endif

	class VRN_CORE_API VortexListPort : public GenericPort<std::vector<Vortex>>
	{
	public:
		VortexListPort( PortDirection direction, const std::string& id, const std::string& guiName = {}, bool allowMultipleConnections = false, Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT );

		virtual Port* create( PortDirection direction, const std::string& id, const std::string& guiName = {} ) const;
		virtual std::string getClassName() const;
		virtual std::string getContentDescriptionHTML() const;
	};
}

#endif // VRN_VORTEXPORT_H
