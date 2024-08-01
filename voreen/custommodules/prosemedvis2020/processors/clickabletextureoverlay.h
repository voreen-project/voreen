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

#ifndef VRN_CLICKABLETEXTUREOVERLAY_H
#define VRN_CLICKABLETEXTUREOVERLAY_H
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <utility>
//header of base class
#include "voreen/core/processors/renderprocessor.h"
//port headers
#include "voreen/core/ports/renderport.h"
#include "../modules/plotting/ports/plotport.h"
//property headers
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/fontproperty.h"
#include "../modules/plotting/properties/colormapproperty.h"
//header of the used shader object
#include "tgt/shadermanager.h"
// CSV-Loader
#include "../utils/samplepointconfigloader.h"
// Stuff
#include "../modules/plotting/datastructures/plotbase.h"



namespace voreen {
/**
 * <Processor Description>
 */
class VRN_CORE_API ClickableTextureOverlay : public RenderProcessor {
public:

    ClickableTextureOverlay();
    virtual Processor* create() const { return new ClickableTextureOverlay(); }
    virtual std::string getClassName() const { return "ClickableTextureOverlay"; }
    virtual std::string getCategory() const { return "Image Processing"; }
    virtual bool isReady() const;

protected:

    virtual void setDescriptions() {
        setDescription("Sample render processor for" \
            "gray-scaling an input image.");
    }
    virtual void process();
    virtual void initialize();
    virtual void deinitialize();
    virtual void onEvent(tgt::Event* e);

private:

    RenderPort inport_;             ///< input of the image which should be modified
    RenderPort outport_;            ///< output of the modified image
    RenderPort concentrationOutport_;
    RenderPort pickingBuffer_;		///< private port used for picking
    PlotPort plotDataInport_;		///< The Port used to supply data for colouring the overlay
    FileDialogProperty fileProp_;
    StringProperty selRegionProp_;	// Property containing the name of the selected region.
    IntProperty selTimestepProp_;		// Number of the selected timestep
    OptionProperty<std::string> selTracerProp_;		// The selected tracer
    ColorMapProperty colorMapProp_;			// ColorMap used in the overlay
    FontProperty heatMapFontProp_;			// Font used in the heatmap
    StringProperty selectedTimeProp_;		// Used to display the time-value of the selected timestep
    StringProperty timeUnitProp_;			// Used to specify the time-unit
    StringProperty measurementUnitProp_;	// Used to specifiy the unit of the measurements
    BoolProperty manualUpperBoundProp_;			// Used to specifiy the upper-bound-mode of the heatmap
    FloatProperty upperBoundProp_;			// Upper Bound for the Heatmap (only used when manualUpperBoundProp_ is set)


    void directoryChanged();
    void mouseEvent(tgt::MouseEvent* e);
    tgt::Shader* shader_;  ///< GLSL shader object used in process()
    std::string mouseDirectory;
    std::vector<SamplePointConfig> regions_;	// Vector containing the regions

    // Map used to buffer the mean-values used in the heatmap.
    // First key is the tracername. Second key is the region name.
    // The vector contains the ordered timesteps.
    std::map<std::string, std::map<std::string, std::vector<plot_t>>> meanBuffer_;
    plot_t maxMeanBuffer_; // Used to buffer the maximum mean-value. (set whenever the plotDataInport_ changes)
    plot_t minMeanBuffer_; // Used to buffer the minmum mean-value. (set whenever the plotDataInport_ changes)

    // Map used to buffer the std-values used in the heatmap.
    // First key is the tracername. Second key is the region name.
    // The vector contains the ordered timesteps.
    std::map<std::string, std::map<std::string, std::vector<plot_t>>> stdBuffer_;
    plot_t maxStdBuffer_; // Used to buffer the maximum std-value. (set whenever the plotDataInport_ changes)
    plot_t minStdBuffer_; // Used to buffer the minimum std-value. (set whenever the plotDataInport_ changes)

    /**
     *	Returns a vector containing the names of all tracers using the column names
     *	from the plotDataInport_
     **/
    std::vector<std::string> getTracerNames() const;

    /**
     *	Returns the mean-value for the given timestep, tracer and region
     *	using the data from plotDataInport_
     *  (first component of the pair is mean, second std)
     **/
    std::pair<plot_t, plot_t> getMeanValue(const int timestepNr, const std::string& tracerName, const std::string& regionName) const;

    /**
     *	Returns the number of timesteps based on the data from plotDataInport_.
     *  In case of an error -1 is returned.
     **/
    int getNumTimesteps() const;

    void renderLegend(const ColorMap& cMap, const float width, const float height, const float xOffset, const float yOffset,
        const plot_t minValue, const plot_t maxValue, const std::string& label);

    void plotDataChanged();

    void drawConcentration();

};

} // namespace


#endif //VRN_CLICKABLETEXTUREOVERLAY_H
