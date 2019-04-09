/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_DYNAMICPYTHONPROCESSOR_H
#define VRN_DYNAMICPYTHONPROCESSOR_H

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/ports/volumeport.h"
#include "../properties/pythonproperty.h"

#include "custommodules/bigdataimageprocessing/properties/interactivelistproperty.h"

namespace voreen {

class DynamicPythonProcessor : public RenderProcessor {
public:
    DynamicPythonProcessor();
    ~DynamicPythonProcessor();

    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "DynamicPythonProcessor"; }
    virtual std::string getCategory() const   { return "Python";                 }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;  }

    virtual bool isReady() const;

    PythonProperty* getPythonProperty() {
        return &pythonProperty_;
    }

protected:

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:

    void onPortListChange();
    void onScriptChange();

    void addPortItem(Port* port);

    std::vector<std::unique_ptr<Port>> portItems_;
    std::map<std::string, std::vector<Port*>> portInstances_;

    InteractiveListProperty portList_;
    BoolProperty enabled_;
    PythonProperty pythonProperty_;
    PythonScript pythonScript_;

    bool valid_;

    static const std::string loggerCat_;
};

} // namespace

#endif // VRN_DYNAMICGLSLPROCESSOR_H
