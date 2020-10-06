/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "parallelcoordinatesselectionproperty.h"

namespace voreen {

ParallelCoordinatesSelectionProperty::ParallelCoordinatesSelectionProperty( const std::string& id, const std::string& guiText, std::vector<int> value, int invalidationLevel, Property::LevelOfDetail lod ) : TemplateProperty( id, guiText, std::move( value ), invalidationLevel, lod )
{}
ParallelCoordinatesSelectionProperty::ParallelCoordinatesSelectionProperty()
    : TemplateProperty()
{}

Property* ParallelCoordinatesSelectionProperty::create() const {
    return new ParallelCoordinatesSelectionProperty();
}
std::string ParallelCoordinatesSelectionProperty::getClassName() const {
    return "ParallelCoordinatesSelectionProperty";
}
std::string ParallelCoordinatesSelectionProperty::getTypeDescription() const {
    return "Vector of Boolean";
}

void ParallelCoordinatesSelectionProperty::serialize( Serializer& serializer ) const {
    Property::serialize( serializer );
}
void ParallelCoordinatesSelectionProperty::deserialize( Deserializer& deserializer ) {
    Property::deserialize( deserializer );
}

}