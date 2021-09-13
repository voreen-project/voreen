#ifndef VRN_LABELVOLUMECREATOR_H
#define VRN_LABELVOLUMECREATOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"

#include "../ports/concretevesselgraphport.h"

namespace voreen {

class VRN_CORE_API LabelVolumeCreator : public Processor {
public:
    LabelVolumeCreator();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "LabelVolumeCreator"; }
    virtual std::string getCategory() const { return "Volume Processing"; }

private:
    VolumePort maskInport_;
    ConcreteVesselGraphPort concreteVesselGraphInport_;

    VolumePort labelVolumeOutport_;

    virtual void process();
};
}
#endif