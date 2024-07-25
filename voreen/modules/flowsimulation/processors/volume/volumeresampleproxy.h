/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMERESAMPLEPROXY_H
#define VRN_VOLUMERESAMPLEPROXY_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/numericproperty.h"

#include "voreen/core/ports/volumeport.h"

namespace voreen {

class VRN_CORE_API VolumeResampleProxy : public Processor,
                                         public PortObserver,
                                         public DataInvalidationObserver {
public:
    VolumeResampleProxy();
    virtual ~VolumeResampleProxy();

    virtual Processor* create() const;

    virtual std::string getClassName() const { return "VolumeResampleProxy"; }

    virtual std::string getCategory() const { return "ProxyVolume"; }

    virtual CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }

    //! @see PortObserver
    virtual void afterConnectionAdded(const Port* source, const Port* connectedPort);
    virtual void beforeConnectionRemoved(const Port* source, const Port*);
    virtual void dataWillChange(const Port* source);
    virtual void dataHasChanged(const Port* source);

    //! @see DataInvalidationObserver
    virtual void dataAboutToInvalidate(const DataInvalidationObservable* data);

protected:
    virtual void setDescriptions() {
        setDescription("This is a resample proxy volume");
    }

    virtual void process();
    void adjustPropertiesToInput();

private:

    VolumePort inport_;
    VolumePort outport_;

    IntVec3Property outputDimensions_;
    ButtonProperty resetResolution_;
    BoolProperty autoResetResolution_;

    static const std::string loggerCat_; ///< category used in logging
};


} // namespace voreen

#endif // VRN_VOLUMERESAMPLEPROXY_H
