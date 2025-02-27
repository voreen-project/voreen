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

#include "parallelcoordinatessource.h"

namespace voreen {

ParallelCoordinatesSource::ParallelCoordinatesSource()
    : Processor()
    , _outport( Port::OUTPORT, "outport", "Parallel Coordinates Axes" )
    , _propertyFileDialog( "property_file_dialog", "File Input", "Select File...", "", "Voreen Parallel Coordinates (*.vpc)", FileDialogProperty::OPEN_FILE )
    , _propertyLoadButton( "property_load_button", "Load" )
{
    addPort( _outport );

    addProperty( _propertyFileDialog );
    addProperty( _propertyLoadButton );
}

Processor* ParallelCoordinatesSource::create() const {
    return new ParallelCoordinatesSource();
}

void ParallelCoordinatesSource::process() {
    if( !_propertyFileDialog.get().empty() ) {
        std::unique_ptr<ParallelCoordinatesAxes> pc(new ParallelCoordinatesAxes(_propertyFileDialog.get()));
        if (pc->members() == 0 || pc->fields() == 0) {
            _propertyFileDialog.set("");
        } else {
            _outport.setData(pc.release(), true);
        }
    }
}

}