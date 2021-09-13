/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2012 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "cmparticleport.h"

namespace voreen {

CMParticlePort::CMParticlePort(PortDirection direction, const std::string& name,
        bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel)
 : GenericPort<CMParticleData>(direction, name, name, allowMultipleConnections, invalidationLevel)
{
    if (isOutport())
        setData(0);
}

CMParticlePort::~CMParticlePort() {
    if ((portData_) && (isOutport())) {
        delete portData_;
        portData_ = 0;
    }
}

tgt::col3 CMParticlePort::getColorHint() const {
	return tgt::col3(18, 183, 92);
}

} // namespace
