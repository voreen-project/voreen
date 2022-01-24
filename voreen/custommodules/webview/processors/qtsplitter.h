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

#ifndef VRN_QTSPLITTER_H
#define VRN_QTSPLITTER_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/network/networkevaluator.h"

#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "modules/base/properties/interactivelistproperty.h"

namespace voreen {

class QtSplitter : public Processor, public ProcessorNetworkObserver {
public:

    enum Orientation {
        HORIZONTAL,
        VERTICAL,
    };

    QtSplitter();
    ~QtSplitter();

    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "QtSplitter";             }
    virtual std::string getCategory() const   { return "WebView";                }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;  }

    virtual bool isReady() const;

    Orientation getOrientation() const;
    std::vector<std::string> getInstances() const;

protected:

    virtual void setDescriptions() {
        setDescription("This processor allows to merge processor widgets similar to the Voreen splitter.");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    virtual void networkChanged();

private:

    OptionProperty<Orientation> orientation_;

    InteractiveListProperty widgets_;
    ButtonProperty updateWidgets_;

    static const std::string loggerCat_;
};

} // namespace

#endif // VRN_QTSPLITTER_H
