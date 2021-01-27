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

#ifndef VRN_TOUCHTABLEBODYPARTSWIDGET_H
#define VRN_TOUCHTABLEBODYPARTSWIDGET_H

#include "modules/touchtable/ports/bodyparts3d/bodypartstextureport.h"

#include "voreen/core/properties/fontproperty.h"

#include "touchtablesnapshotwidget.h"
#include "touchtablemenuwidgetscrollable.h"
#include "../touchtablecontrolelement.h"
#include "../touchtablepreviewicon.h"
#include "../touchtableoverlay.h"

namespace voreen{

    typedef SnapshotPreview BodyPartPreview;

    class TouchTableBodyPartsWidget : public TouchTableMenuWidgetScrollable {
    public:
        TouchTableBodyPartsWidget();

        ~TouchTableBodyPartsWidget();

        Processor* create() const{
            return new TouchTableBodyPartsWidget();
        }

        virtual bool isReady() const{
            return port_.isConnected();
        }

        virtual std::string getClassName() const{
            return "TouchTableBodyPartsWidget";
        }

        virtual void handleTouchPoints(const std::deque<tgt::TouchPoint>& tp);

    protected:

        void setDescriptions(){
            setDescription("manages removed BodyParts");
        }

        virtual void renderComponents();

        virtual void initialize();

        virtual void deinitialize();

        virtual void updateComponents();

        virtual void updateBodyPartsMenu();

        virtual void updateScrollableMenuPosition();

        virtual void handleScrollableMenuSelection(std::string filename);

        virtual void process();

        bool state_;
        int selected_;
        float sliderPos_;

        tgt::Texture* restoreTex_;
        tgt::Texture* defaultTex_;

        BodyPartPreview preview_;

        TouchTableControlElement restoreElem_;

        IntProperty bodypartsMenuHeight_;
        IntProperty bodypartsMenuWidth_;
        IntProperty restoreSignal_;

        std::vector<TouchTablePreviewIcon> bodyParts_;

        TouchTableScrollableMenu bodyPartsMenu_;

        BodyPartsTexturePort inPort_;

        static const std::string loggerCat_;
    };
};
#endif //VRN_TOUCHTABLEBODYPARTSWIDGET_H
