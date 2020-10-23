

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/geometryport.h"

#include "modules/flowanalysis/ports/vortexport.h"


namespace voreen
{
	class VortexSelector : public Processor
	{
	public:
		VortexSelector();
		virtual Processor *create() const { return new VortexSelector(); }
		virtual std::string getClassName() const { return "VortexSelector"; }
		virtual std::string getCategory() const { return "Vortex Processing"; }
		bool isReady() const override
		{
			return _inVortexList.isReady();
		}

	protected:
		virtual void process();

	private:
		VortexListPort _inVortexList;
		VortexPort _outVortex;
		GeometryPort _outCoreline;
		IntProperty _selectedIndex;
	};
} // namespace voreen
