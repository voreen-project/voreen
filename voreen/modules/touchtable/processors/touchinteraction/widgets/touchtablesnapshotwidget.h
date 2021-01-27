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

#ifndef VRN_TOUCHTABLESNAPSHOTWIDGET_H
#define VRN_TOUCHTABLESNAPSHOTWIDGET_H

#include "../touchtableoverlay.h"
#include "touchtablemenuwidget.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/fontproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/filedialogproperty.h"

#include "../touchtablecontrolelement.h"
#include "../touchtablemenuframe.h"
#include "../touchtablepreviewicon.h"
#include "touchtablemenuwidgetscrollable.h"

#include "tgt/filesystem.h"

#include "../touchtableslider.h"

namespace voreen {

    /**
    * A struct to encapsulate position, and texture of the snapshot Preview
    *
    *@member posLL_ lower left position of the preview
    *@member posUR_ upper right position of the preview
    *@member previewTexture_ texture of the preview to be rendered
    */
    struct SnapshotPreview{

        SnapshotPreview(tgt::ivec2 posLL, tgt::Texture* previewTexture)
            : posLL_(posLL)
            , previewTexture_(previewTexture)
        {
        }

        tgt::ivec2 posLL_;
        tgt::ivec2 posUR_;
        tgt::Texture* previewTexture_;
    };

    /**
    * A struct to encapsulate the message settings
    *
    *@memb message_ message to be displayed
    *@memb display_ true if message is displayed, false otherwise
    */
      struct Message {

            Message() {
                message_ = "";
                display_ = false;
            }

            std::string message_;
            bool display_;
        };

    /**
    * TouchTableSnapShotWidget inherits from TouchTableMenuWidget. When active through this widget
    * screenshots can be taken and managed. Taking the screenshots and saving them in the file system is
    * done through the respective functionalities of the canvas
    */
    class TouchTableSnapShotWidget : public TouchTableMenuWidgetScrollable {

    public:

        TouchTableSnapShotWidget();

        ~TouchTableSnapShotWidget();

        Processor* create() const {
            return new TouchTableSnapShotWidget();
        }

        virtual bool isReady() const {
            return port_.isConnected();
        }

        virtual std::string getClassName() const {
            return "TouchTableSnapShotWidget";
        }

        /**
         * Handles touch points within the snap shot overlay menu
         *
         * @param tp a std::deque containing all the touch points that should be handled by the widget.
         * @param overlay the TouchTableOverlay processor which has been catching the touch point events before.
         */
        virtual void handleTouchPoints(const std::deque<tgt::TouchPoint>& tp);

        /** To handle invalid snapshot directory */
        virtual void pressed() override;

    protected:

        void setDescriptions(){
            setDescription( "A TouchTableModeWidget through which Snapshots of the current rendering can be taken and managed in a preview list");
        }

        /**
         * Renders the snap shot overlay menu, the preview of the screenshots, control elements etc.
         */
        virtual void renderComponents();

        virtual void initialize();

        virtual void deinitialize();

         /**
         * Is called if properties change to update the control element attributes.
         */
        virtual void updateComponents();

        /**
        * updates the content of the scrollable menu with its previews of the taken snapshots
        */
        virtual void updateSnapshotMenu();
        /**
        *updates the coordintes of the scrollable menu with its previews of the taken snapshots (overwritten from base class)
        */
        virtual void updateScrollableMenuPosition();

        /**
        *shows the currently selected snapshot in the bigger preview
        *@param filename filename of the selected snapshot
        */
        virtual void handleScrollableMenuSelection(std::string filename);

        /**
        * sets the message, showing while loading
        *@param dialog dialog that is being shown while snapshot is being taken
        *@param show true if dialog is supposed to be shown, false otherwise
        */
        virtual void setLoadingMessage(std::string dialog, bool show);

        /**
        * saves the taken snap shot and loads the texture into the snapshot preview
        *@param filename filename of the taken screenshot with file ending
        */
        virtual void saveSnapshot(std::string filename);

         /**
         * Handles timer events for selecting default camera orientations
         *@ te given timer event
         */
        virtual void timerEvent(tgt::TimeEvent* te);

        //textures for control elements
        tgt::Texture* snapshotTex_, * deleteTex_;

        TouchTableControlElement snapshotElem_;                ///< control element to take snapshots
        TouchTableControlElement deleteElem_;                ///< control element to delete snapshots

        //properties
        FileDialogProperty snapshotDirectory_;                ///< saves the directory path, where screenshots are saved
        IntProperty snapshotMenuHeight_;                    ///< Height of the menu containing small icons of all taken snapshots
        IntProperty snapshotMenuWidth_;                        ///< Length of the menu containing small icons of all taken snapshots
        IntProperty snapshotIconSize_;                        ///< size of small icons representing the taken screenshots
        FontProperty fontProperty_;                            ///< font options for the gallery of snapshot previews


        std::vector<TouchTablePreviewIcon> snapshots_;        ///< vector saving all previewicons of all screenshots being taken

        SnapshotPreview snapshotPreview_;                    ///< bigger preview of selecte snapshot from the scrollable menu

        std::string currentFilename_;                        ///< filename of the currently taken screenshot

        Message fileMessage_;                                ///< message to be displayed while saving/loading

        tgt::Timer* messageTimer_;                            ///< timer object for displaying save/load messages

        tgt::Timer* saveSnapshotTimer_;                        ///< time object for saving the taken screenshot

        tgt::EventHandler eventHandler_;                    ///< a local EventHandler that is added to the timer


        static const std::string loggerCat_;
    };

} // namespace

#endif // VRN_TOUCHTABLESNAPSHOTWIDGET_H

