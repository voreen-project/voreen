#include "vortexcollectionsource.h"

namespace voreen
{
	VortexCollectionSource::VortexCollectionSource() : Processor(),
		_outportVortexCollection( Port::OUTPORT, "outport", "Vortex Collection" ),
		_propertyFileDialog( "property_file_dialog", "File Input", "Select File...", "", "Vortex Collection (*.vc)", FileDialogProperty::OPEN_FILE, Processor::VALID ),
		_propertyLoadButton( "property_load_button", "Load" )
	{
		this->addPort( _outportVortexCollection );

		this->addProperty( _propertyFileDialog );
		this->addProperty( _propertyLoadButton );
	}

	Processor* VortexCollectionSource::create() const
	{
		return new VortexCollectionSource();
	}

	std::string VortexCollectionSource::getClassName() const
	{
		return "VortexCollectionSource";
	}
	std::string VortexCollectionSource::getCategory() const
	{
		return "Vortex Extraction";
	}

	void VortexCollectionSource::process()
	{
		if( _propertyFileDialog.get() != "" )
		{
			if( auto stream = std::ifstream( _propertyFileDialog.get(), std::ios::in | std::ios::binary ) )
				_outportVortexCollection.setData( new VortexCollection( stream ) );
		}
	}
}
