#ifndef VRN_PROSEMEDVIS2020MODULE_H
#define VRN_PROSEMEDVIS2020MODULE_H

#include "voreen/core/voreenmodule.h"

namespace voreen {

class ProseMedvis2020Module : public VoreenModule {
public:
    ProseMedvis2020Module(const std::string& modulePath);

    virtual std::string getDescription() const {
        return "Module for WS2021 ProSem";
    }
};

}


#endif
