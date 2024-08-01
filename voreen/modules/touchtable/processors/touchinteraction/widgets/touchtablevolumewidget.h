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

#ifndef VRN_TOUCHTABLEVOLUMEWIDGET_H
#define VRN_TOUCHTABLEVOLUMEWIDGET_H

#include "../touchtableoverlay.h"
#include "touchtablemenuwidget.h"
#include "touchtablemenuwidgetscrollable.h"

#include "voreen/core/properties/fontproperty.h"
#include "voreen/core/properties/volumeurlproperty.h"
#include "voreen/core/properties/volumeinfoproperty.h"

#include "tgt/filesystem.h"

namespace voreen{

    class TouchTableVolumeWidget : public TouchTableMenuWidgetScrollable {

    public:
        TouchTableVolumeWidget();

        Processor* create() const {
            return new TouchTableVolumeWidget();
        }

        virtual std::string getClassName() const {
            return "TouchTableVolumeWidget";
        }

        /** Deactivtes, if the volume directory is set invalid. */
        virtual void pressed() override;

        /**
         * Handles touch points within the volume menu
         *
         * @param tp a std::deque containing all the touch points that should be handled by the widget.
         */
        virtual void handleTouchPoints(const std::deque<tgt::TouchPoint>& tp);

        /**
         * sets linked volumeSource with url build from volume directory and file string
         *
         * @param file selection of scrollablemenu (name of volume)
         */
        virtual void handleScrollableMenuSelection(std::string file);

        virtual void setValidFileExtensions(std::vector<std::string> extensions);
        virtual std::vector<std::string> getValidFileExtensions();

    protected:

        void setDescriptions(){
            setDescription("loads volume data");
        }

        virtual void initialize();

        virtual void deinitialize();

        virtual void renderComponents();
        virtual void updateComponents();
        virtual void updateScrollableMenuPosition();
        virtual void updateScrollableMenuContent();
        virtual void updateFontProp();

        /**
        * returns true if an extension is valid (is in validFileExtension vector)
        *
        * @param extension the extension to check
        */
        virtual bool isExtensionValid(std::string extension);

        std::vector<std::string> validFileExtensions_;
        FileDialogProperty volumesDirectory_;               ///< property for selecting a directory in which to search for TF presets
        FontProperty fontProp_;                             ///< property for font rendering of preset menu
        VolumeURLProperty volumeURL_;                        ///< property to set volume source
        VolumeInfoProperty volumeInfo_;                        ///< property linked with voulmeURL

        int elementsToRender_;                                ///< number of elements rendered into content menu
        static const std::string loggerCat_;
    };

} // namespace

#endif
