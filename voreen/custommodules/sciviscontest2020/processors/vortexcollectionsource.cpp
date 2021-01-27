/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "vortexcollectionsource.h"

namespace voreen {
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
