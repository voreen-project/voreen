/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_COMPONENTIDSELECTOR_H
#define VRN_COMPONENTIDSELECTOR_H

#include "voreen/core/processors/imageprocessorbypassable.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/properties/string/stringlistproperty.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "voreen/core/ports/volumeport.h"

#include "tgt/font.h"
#include "tgt/glmath.h"
#include "tgt/immediatemode/immediatemode.h"

namespace voreen {

    class VRN_CORE_API ComponentIdSelector : public ImageProcessorBypassable {
    public:
        ComponentIdSelector();
        ~ComponentIdSelector();
        virtual Processor* create() const;

        virtual std::string getCategory() const  { return "Utility";         }
        virtual std::string getClassName() const { return "ComponentIdSelector"; }
        virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL; }
        virtual bool isUtility() const           { return true; }

        virtual bool isReady() const;

        void serialize(Serializer& s) const override;
        void deserialize(Deserializer& s) override;

    protected:
        virtual void setDescriptions() {
            setDescription(
                    "Lorem Ipsum"
            );
        }

        void process();
        virtual void onEvent(tgt::Event* e);

    private:

        void fillList();

        RenderPort raycastingResultInport_;
        RenderPort fhpInport_;
        VolumePort refInport_;
        RenderPort outport_;

        IntProperty maxID_;
        IntProperty selectedID_;
        StringListProperty components_;
        IntVec3Property selectedPosition_;
        IntProperty idSearchRadius_;

        bool mouseDown_;
        tgt::ivec2 mousePosition2D_;

        std::set<int> selection_;

        void mouseEvent(tgt::MouseEvent* e);
        tgt::ivec2 clampToViewport(tgt::ivec2 mousePos);
    };

} // namespace

#endif // VRN_COMPONENTIDSELECTOR_H
