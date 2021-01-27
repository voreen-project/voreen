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

#ifndef VRN_LINEPROFILE_H
#define VRN_LINEPROFILE_H
//general defines
#include "voreen/core/processors/imageprocessor.h"
//properties
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/optionproperty.h"
//ports
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/textport.h"
#include "../../plotting/ports/plotport.h"
//plotting
#include "../../plotting/datastructures/plotdata.h"
#include "../../properties/plotselectionproperty.h"
#include "../../properties/plotentitiesproperty.h"
//for rotation matrix
#include "tgt/glmath.h"
//utils
#include "../utils/LineProfile/lp_graphic.h"
#include "../utils/LineProfile/lp_measure.h"
#include "../utils/LineProfile/lp_texteditor.h"
#include "../utils/LineProfile/lp_plotmodifier.h"
#include "../utils/LineProfile/lp_maxthreshold.h"
#include "../utils/LineProfile/lp_eigenfunctionfit.h"
//textoutput
#include <sstream>
#include <vector>

//setting namespace
namespace voreen {
 /*
  * class for drawing a line into a volume slice and measure the intensity along this line
  */
class VRN_CORE_API LineProfile : public ImageProcessor {
public:
    LineProfile();
    ~LineProfile();
    virtual Processor* create() const;

    virtual std::string getCategory() const  { return "Utility"; }
    virtual std::string getClassName() const { return "LineProfile"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL; }
    virtual bool isUtility() const           { return false; }

    virtual bool isReady() const;
    virtual void invalidate(int inv = 1);

    //The event functions
    void drawLine(tgt::MouseEvent* e);        //for drawing the line
    void moveSphere(tgt::MouseEvent* e);    //for moving the spheres

protected:
    virtual void setDescriptions() {
        setDescription("");
    }

    virtual void interactionModeToggled();
    virtual void process();
    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);

private:
    //if line is drawn
    bool hasLine_;
//ports
    VolumePort volSingleInport_;//must be connected, function fit on this volume
    VolumePort volMultiInport_;    //optional volumes shown in plot

    RenderPort imgInport_;
    RenderPort fhpInport_;
    RenderPort imgOutport_;
    RenderPort imgTempport_; //private port

    PlotPort plotOutport_;
    TextPort textOutport_;

//mouse events
    EventProperty<LineProfile>* mouseEventDrawLine_;    //event to draw the arrow
    EventProperty<LineProfile>* mouseEventMoveSphere_;    //event for showing sphere by move over
    EventProperty<LineProfile>* mouseEventPressSphere_; //event to grap sphere
    bool mouseDownDrawLine_;                            //if mouse is pressed for drawing the line
    bool mouseDownMoveSphere_;                            //if mouse is pressed for moving a sphere
    enum sphere_enum {NON, STA, MID, END} graped_;        //shows, if a sphere is graped
    //crop mouse position to canvas size
    tgt::ivec2 cropToViewport(tgt::ivec2 mousePos);

//general properties and settings
    enum calculateNew {
        CN_NON,
        CN_PLOTDATA,
        CN_FITFUNCTION,
        CN_THRESHOLD,
        CN_TEXTOUTPUT,
        CN_IMAGE
        } calculateNew_;
    CameraProperty camera_;
    void preSettings();
    //sets hasChanged_ = ..., if a property has changed
    void plotdataMakeNew();
    void fitFunctionMakeNew();
    void thresholdMakeNew();
    void textOutputMakeNew();
    void imageMakeNew();
    bool inportVolume(std::vector<const VolumeBase* >& vvh); //inport volume from inport

//plotting
    PlotData* pVolData_;                        //plotdata form input
    PlotData* pModData_;                        //plotdata modified
    OptionProperty<LP_Measure::fetchOption> fetchOption_;    //which voxel fetch is being used
    FloatProperty samplingRate_;                //the sampling rate
    FloatProperty interactionQuality_;            //the factor to modify the sampling rate
    bool setPlotData(std::vector<const VolumeBase* >& vvh, PlotData& pData);

//lm_algorithm
    BoolProperty calculateLM_;  //true, if function fit is used
        void calculateLMOnChange();
    BoolProperty useLM_;        //true, if threshold uses ff
    BoolProperty singleRepresentation_;
    OptionProperty<LP_EigenFunctionFit::fitFunctionMode> fitOption_;
        void fitOptionOnChange();
    IntProperty fitOptionAddOn_;//additional settings for fitOption_
    StringProperty fitRepresentation_;
    //calculate Fit unction
    void setFitFunction();

//maxthreshold
    float maxValue_;            //active value of maximum
    float threshold_;            //active value of tLine
    std::vector< std::pair< float,float > > thresholdArea_;        //threshold visualisation
    //sets all members
    void setMaxValue();
    void setThresholdValue();
    void setThresholdLine();

    BoolProperty calculateThreshold_;        //true, if threshold calculation is used
        void calculateThresholdOnChange();
    BoolProperty showIntervals_;    //show intervals
    BoolProperty useAbsLine_;        //use absolute or percent values
        void useAbsLineOnChange();
    FloatProperty thresholdLineAbs_;    //the absolute value
        void setLineAbsMax();            //set max value of property
    FloatProperty thresholdLinePer_;    //the percent value

//rendering
    //arrow
    BoolProperty arrowVisible_;                //if arrow should be rendered
    FloatProperty arrowSize_;                //the size/width of the arrow
    FloatVec4Property arrowColor_;            //the color of the arrow
    FloatVec4Property arrowMarkColor_;        //the color the marked sections
    FloatVec4Property arrowHotColor_;        //the color of sections above a threshold
    PlotSelectionProperty selectionProp_;    //linked property for marked visualisation
    //draw the arrow
    void drawLine();

    //move spheres
    BoolProperty sphereMoveVisible_;        //if move spheres are always visible
    FloatProperty sphereMoveSize_;            //the size of the move spheres
    FloatVec4Property sphereMoveColor_;        //the color of the move spheres
    //Positions
        tgt::vec3 staSpherePos_; tgt::vec3 midSpherePos_; tgt::vec3 endSpherePos_;
        //true, if the mouse is over the sphere of visible = true
        bool drawStaSphere_; bool drawMidSphere_; bool drawEndSphere_;

    float sphereMaxPosition_;    //the position of the max sphere
    BoolProperty showMax_;                    //show the maximum value of the plotdata
    FloatProperty sphereMaxSize_;            //the size of the max sphere
    FloatVec4Property sphereMaxColor_;        //the color of the max sphere


//textoutput
    std::stringstream strstr_; //used for textoutput
    void setTextOutport();

};

} //end namespace

#endif //end VRN_LineProfile_H
