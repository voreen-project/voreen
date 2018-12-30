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

#ifndef VRN_RADARGLYPHRENDERER2D_H
#define VRN_RADARGLYPHRENDERER2D_H

//super class
#include "radarglyphrendererbase.h"

#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/interaction/mwheelnumpropinteractionhandler.h"


namespace voreen {

class RadarGlyphRenderer2D : public RadarGlyphRendererBase
{
public:
    RadarGlyphRenderer2D();
    virtual ~RadarGlyphRenderer2D();

    virtual Processor* create() const { return new RadarGlyphRenderer2D(); }
    virtual std::string getClassName() const { return "RadarGlyphRenderer2D"; }
    virtual void process();

private:
    enum RenderMode2D {
        NORMAL_RADAR_GLYPHS_2D_XY = 0,
        NORMAL_RADAR_GLYPHS_2D_YZ = 1,
        NORMAL_RADAR_GLYPHS_2D_XZ = 2
    };
    friend class OptionProperty<int>;

    IntProperty sliceIndex_;
    ColorProperty boundaryColorProp_;

    //2D Control
    FloatProperty zoomProp_;
    FloatProperty xOffsetProp_;
    FloatProperty yOffsetProp_;
    EventProperty<RadarGlyphRenderer2D>* mouseEventShift_;
    void shiftEvent(tgt::MouseEvent* e);
    MWheelNumPropInteractionHandler<int> mwheelCycleHandler_;
    MWheelNumPropInteractionHandler<float> mwheelZoomHandler_;
    tgt::ivec2 mousePosition_;

    //onChange
    void renderModeOnChange();
    void inportOnChange();
    void sliceIndexOnChange();

    //fill VBOs
    void fillGlyphVBO();
    //void fillIndexVBO();
    tgt::ivec3 permute(int slice, int y, int x);

    //members needed for 2d projection
    size_t x_count_, y_count_; /// number of glyphs in this direction
    float x_offset_, y_offset_;/// pixel offset in this direction
    float glyphMaxDiameter_;  /// pixelsize of the glyph
};

}   // namespace

#endif  // VRN_RADARGLYPHRENDERER2D_H
