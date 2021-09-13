#ifndef VRN_TEMPLATEVESSELGRAPHLOADER_H
#define VRN_TEMPLATEVESSELGRAPHLOADER_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"

#include "../ports/templatevesselgraphport.h"
#include "../datastructures/templatevesselgraph.h"

namespace voreen {

class VRN_CORE_API TemplateVesselGraphLoader : public Processor {
public:
    TemplateVesselGraphLoader();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "TemplateVesselGraphLoader"; }
    virtual std::string getCategory() const { return "Volume Processing"; } //?

private:
    TemplateVesselGraphPort vesselGraphOutput_;

    virtual void process();
};
}

#endif
