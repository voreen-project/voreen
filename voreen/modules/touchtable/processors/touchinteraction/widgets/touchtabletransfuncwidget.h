/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_TOUCHTABLETRANSFUNCWIDGET_H
#define VRN_TOUCHTABLETRANSFUNCWIDGET_H

#include "../touchtableoverlay.h"
#include "touchtablemenuwidgetscrollable.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/fontproperty.h"

#include "voreen/core/datastructures/transfunc/1d/1dkeys/transfunc1dkeys.h"
#include "voreen/core/datastructures/volume/histogram.h"

#include "../touchtablecontrolelement.h"
#include "../touchtablemenuframe.h"

#include "tgt/filesystem.h"

#include "../touchtableslider.h"

namespace voreen {

    /**
    * Manages a TF menu in which the histogram and transfer function are displayed, the TF domain may be adjusted and TF presets may be loaded.
    */
    class TouchTableTransFuncWidget : public TouchTableMenuWidgetScrollable {

    public:

        TouchTableTransFuncWidget();

        Processor* create() const {
            return new TouchTableTransFuncWidget();
        }

        virtual bool isReady() const {
            return port_.isConnected();
        }

        virtual std::string getClassName() const {
            return "TouchTableTransFuncWidget";
        }

        virtual void handleTouchPoints(const std::deque<tgt::TouchPoint>& tp);

        virtual void handleScrollableMenuSelection(std::string file);

    protected:

        void setDescriptions(){
            setDescription("Manages a TF menu in which the histogram and transfer function are displayed, \
                    the TF domain may be adjusted and TF presets may be loaded.");
        }

        /**
         * Renders the control elements inside the transfer function overlay menu as well as the histogram and TF.
         */
        virtual void renderComponents();

        virtual void initialize();

        virtual void deinitialize();

        /**
         * Updates the histogram if the input volume changed and calls process() of parent class.
         * @see TouchTableMenuWidgetScrollable
         */
        virtual void process();

        /**
         * Renders the histogram and transfer function on top of a checkerboard background.
         *
         * @param ll lower left (in canvas, ie. pixel, coordinates)
         * @param ur upper right (in canvas, ie. pixel coordinates)
         * @param viewportSize size of the whole viewport (eg. canvas connected to TouchTableOverlay)
         * @param renderBackground determines if the checkerboard background should be rendered underneath
         * @param renderTF determines if the TF should be rendered on top of the background
         * @param renderHistogram determines if the histogram should be rendered on top of the background and TF
         */
        virtual void renderHistogram(tgt::ivec2 ll, tgt::ivec2 ur, tgt::ivec2 viewportSize, bool renderBackground, bool renderTF, bool renderHistogram);

        virtual void updateComponents();

        virtual void updateScrollableMenuPosition();

        /**
         * The preset menu has either been opened or closed.
         * Updates the status of the corresponding control element and (if necessary) the vector with TF preset file names.
         */
        virtual void updateScrollableMenuContent();

        /// zoom in on transfer function and histogram
        virtual void zoomIn();

        /// zoom out on transfer function and histogram
        virtual void zoomOut();

        /// reset zoom on transfer function and histogram
        virtual void resetZoom();

        virtual void adaptSlidersToZoom();
        virtual void adaptSliderPositionsToDomain();

        TransFunc1DKeysProperty transfunc_;                       ///< transfer function property (to be linked with raycaster and other processors in the voreen network)

        TransFunc1DKeysProperty internalTransFunc_;               ///< transfer function for temporarily loading the color maps opf transfer functions

        VolumePort volumePort_;                             ///< volume port (needed for computing the histogram)

        FileDialogProperty presetDirectory_;                ///< property for selecting a directory in which to search for TF presets

        BoolProperty renderBackground_;                     ///< if enabled the background for the histogram is rendered
        BoolProperty renderTF_;                             ///< if enabled the transfer function is rendered above the background
        BoolProperty renderHistogram_;                      ///< if enabled the histogram is rendered above the background and tf

        FloatVec2Property viewLeftRight_;                   ///< view positions for zooming and painting the histogram

        //control elements that control whhich parts of the histogram are rendered
        TouchTableControlElement backgroundElem_;           ///< used to control the background bool property
        TouchTableControlElement tfElem_;                   ///< used to control the tf bool property
        TouchTableControlElement histogramElem_;            ///< used to control the histogram bool property

        //control elements for zooming on the histogram
        TouchTableControlElement zoomIn_;
        TouchTableControlElement zoomOut_;
        TouchTableControlElement zoomReset_;

        TouchTableControlElement fitToData_;                ///< fit TF to input volume data

        TouchTableControlElement loadPreset_;               ///< control element for loading TF presets

        TouchTableTwoIndicatorSlider domainSlider_;         ///< slider for handling TF domain

        //textures for control elements
        tgt::Texture* gridTex_, * mappingTex_, * histTex_;
        tgt::Texture* zoomInTex_, * zoomOutTex_, * zoomResetTex_;
        tgt::Texture* fitToDataTex_;
        tgt::Texture* presetTex_;

        Histogram1D* histogram_;                            ///< the histogram of the current input volume

        std::vector<std::pair<std::string, TransFunc1DKeys*> > tfPresetFiles_;            ///< contains all tf preset file names availabe in the preset menu

        std::vector<int> domainSliderIDs_;                  ///< IDs of touch points that are associated with the domain slider

        static const std::string loggerCat_;
    };

} // namespace

#endif // VRN_TOUCHTABLETRANSFUNCWIDGET_H
