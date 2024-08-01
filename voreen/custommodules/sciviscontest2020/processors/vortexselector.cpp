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

#include <chrono>

#include "vortexselector.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"

namespace voreen
{
	VortexSelector::VortexSelector() : Processor(),
											   _inVortexList(Port::INPORT, "inVortexList", "Vortex list"),
											   _outVortex(Port::OUTPORT, "outVortex", "Vortex"),
											   _outCoreline(Port::OUTPORT, "outCoreline", "Coreline"),
											   _selectedIndex("selectedIndex", "Index", 0, 0, std::numeric_limits<int>::max())
	{
		this->addPort(_inVortexList);
		this->addPort(_outVortex);
		this->addPort(_outCoreline);
		this->addProperty(_selectedIndex);
		_inVortexList.onNewData(LambdaFunctionCallback([this] {
			_selectedIndex.setMaxValue(std::max(0, static_cast<int>(_inVortexList.getData()->size()) - 1));
		}));
	}

	void VortexSelector::process()
	{
		if (!_inVortexList.hasData() || !_inVortexList.getData()->size())
			return;

		const auto outVortex = new Vortex((*_inVortexList.getData())[_selectedIndex.get()]);
		auto outCoreline = new PointListGeometryVec3();
		outCoreline->setData(outVortex->coreline());

		_outVortex.setData(outVortex);
		_outCoreline.setData(outCoreline);
	}
} // namespace voreen
