#include <chrono>

#include "vortexselector.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

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
