/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_LINKEVALUATOCOLORSWITCH_H
#define VRN_LINKEVALUATOCOLORSWITCH_H

#include "voreen/core/properties/link/linkevaluatorbase.h"
#include "voreen/core/properties/color/colorswitchproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/boolproperty.h"

namespace voreen {


class LinkEvaluatorColorSwitchBool : public LinkEvaluatorBase{
public:
    virtual void eval(Property* src, Property* dst);
    virtual bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
    virtual std::string getClassName()  const { return "LinkEvaluatorColorSwitchBool";     }
    virtual std::string getGuiName()    const { return "Color Switch <-> Bool";            }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorColorSwitchBool(); }
};
/*
class LinkEvaluatorColorSwitchColor : public LinkEvaluatorBase{
public:
    virtual void eval(Property* src, Property* dst);
    virtual bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
    virtual std::string getClassName()  const { return "LinkEvaluatorColorSwitchColor";     }
    virtual std::string getGuiName()    const { return "Color Switch <-> Color";            }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorColorSwitchColor(); }
};*/

} // namespace

#endif // VRN_LINKEVALUATORBOOLINVERT_H
