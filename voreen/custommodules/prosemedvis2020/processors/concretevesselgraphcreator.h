#ifndef VRN_CONCRETEVESSELGRAPHCREATOR_H
#define VRN_CONCRETEVESSELGRAPHCREATOR_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/volumeport.h"
#include "modules/vesselnetworkanalysis/ports/vesselgraphport.h"
#include "../ports/templatevesselgraphport.h"
#include "../ports/concretevesselgraphport.h"

#include "../datastructures/templatevesselgraph.h"
#include "../datastructures/concretevesselgraph.h"

#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
namespace voreen {
    
    class VRN_CORE_API ConcreteVesselGraphCreator : public Processor {
        public:
            ConcreteVesselGraphCreator();
            virtual Processor* create() const;

            virtual std::string getClassName() const { return "ConcreteVesselGraphCreator"; }
            virtual std::string getCategory() const { return "Volume Processing"; }

        private:
            VolumePort flowInport;
            TemplateVesselGraphPort templateVesselGraphInport;
            VesselGraphPort vesselGraphInport;

			FileDialogProperty filePath;
			ButtonProperty saveButton;

            ConcreteVesselGraphPort concreteVesselGraphOutport;

            virtual void process();
			void saveCurrentConcreteVesselGraph();
    };
}

#endif
