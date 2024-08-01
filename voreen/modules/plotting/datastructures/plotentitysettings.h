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

#ifndef VRN_PLOTENTITYSETTINGS_H
#define VRN_PLOTENTITYSETTINGS_H

#include "plotdata.h"
#include "colormap.h"
#include "voreen/core/io/serialization/serializable.h"
#include "tgt/vector.h"

namespace voreen {

/**
 * This serializable class stores information about a plot entity like a line,
 * a group of candle sticks or bars etc.
 */
class VRN_CORE_API PlotEntitySettings : public Serializable {
public:

    // enum for specifying entity
    enum Entity {
        NONE = -1,
        LINE = 0,
        BAR = 1,
        SURFACE = 2,
        SCATTER = 3
    };

    /// enum for the line style
    enum LineStyle {
        CONTINUOUS = 0,
        DOTTED = 1,
        DASHED = 2
    };

    enum PolygonFillStyle {
        FILLNONE = 0,
        NONZERO = 1,
        EVENODD = 2,
        TEXTURE = 3
    };

    /// enum for the glyph style
    enum GlyphStyle {
        POINT = 0,    ///< simple point
        CIRCLE = 1,   ///< circle or sphere
        TRIANGLE = 2, ///< triangle or tetraeder
        QUAD = 3      ///< rectangle or cube
    };

    PlotEntitySettings();

    /// constructor for line
    PlotEntitySettings(int lineCI, tgt::Color lineColor, LineStyle style,
        bool splineFlag, int errorCI, tgt::Color errorColor,  bool errorbarFlag);

    /// constructor for bar
    PlotEntitySettings(int barsCI, tgt::Color barColor);

    /// constructor for surface
    PlotEntitySettings(int surfaceCI, tgt::Color surfaceColor, bool wireFlag, tgt::Color wireColor,
        ColorMap cm, bool heigtmapFlag, int optionalCI);

    /// constructor for scatter
    PlotEntitySettings(int scatterCI, tgt::Color color, float glyphSizeMin ,float glyphSizeMax, GlyphStyle style);

    /**
    * @see Property::serialize
    */
    virtual void serialize(Serializer& s) const;

    /**
    * @see Property::deserialize
    */
    virtual void deserialize(Deserializer& s);

    /// return true, whether column index and column label match
    bool fitsTo(const PlotData* data) const;

    Entity getEntity() const;
    void setEntity(Entity entity);

    int getMainColumnIndex() const;
    void setMainColumnIndex(int value);

    int getCandleTopColumnIndex() const;
    void setCandleTopColumnIndex(int value);

    int getCandleBottomColumnIndex() const;
    void setCandleBottomColumnIndex(int value);

    int getStickTopColumnIndex() const;
    void setStickTopColumnIndex(int value);

    int getStickBottomColumnIndex() const;
    void setStickBottomColumnIndex(int value);

    int getOptionalColumnIndex() const;
    void setOptionalColumnIndex(int value);

    int getSecondOptionalColumnIndex() const;
    void setSecondOptionalColumnIndex(int value);

    bool getCandleStickFlag() const;
    void setCandleStickFlag(bool value);

    LineStyle getLineStyle() const;
    void setLineStyle(LineStyle value);

    bool getSplineFlag() const;
    void setSplineFlag(bool value);

    bool getErrorbarFlag() const;
    void setErrorbarFlag(bool value);

    bool getWireOnlyFlag() const;
    void setWireOnlyFlag(bool value);

    bool getHeightmapFlag() const;
    void setHeightmapFlag(bool value);

    ColorMap getColorMap() const;
    void setColorMap(ColorMap colorMap);

    float getMinGlyphSize() const;
    void setMinGlyphSize(float value);

    float getMaxGlyphSize() const;
    void setMaxGlyphSize(float value);

    GlyphStyle getGlyphStyle() const;
    void setGlyphStyle(GlyphStyle value);

    tgt::Color getFirstColor() const;
    void setFirstColor(tgt::Color value);

    tgt::Color getSecondColor() const;
    void setSecondColor(tgt::Color value);

    std::string getTexturePath() const;
    void setTexturePath(std::string path);

    bool getUseTextureFlag() const;
    void setUseTextureFlag(bool value);

    bool operator==(const PlotEntitySettings& rhs) const;

private:

    /// type of the entity
    Entity entity_;

    // an entity stores one or more column indices and labels
    int mainColumnIndex_;           ///< column index for main axis
    int candleTopColumnIndex_;      ///< column index for top of candle
    int candleBottomColumnIndex_;   ///< column index for bottom of candle
    int stickTopColumnIndex_;       ///< column index for top of candle stick
    int stickBottomColumnIndex_;    ///< column index for bottom if candle stick
    int optionalColumnIndex_;       ///< column index of optional column, the value -1 means, that the column is explicitly -not- selected
    int secondOptionalColumnIndex_; ///< another optional column index

    tgt::Color firstColor_;         ///< main color of the entity
    tgt::Color secondColor_;        ///< secondary color of the entity
    ColorMap colorMap_;             ///< colormap to apply

    // line plot related settings:
    LineStyle lineStyle_;           ///< drawing style of line
    bool splineFlag_;               ///< flag whether this is a spline
    bool errorbarFlag_;             ///< flag whether the error is drawn as error bars
    bool candleStickFlag_;          ///< flag whether this is a candle stick

    // surface plot related settings:
    bool wireOnlyFlag_;             ///< flag whether only the wireframe will be rendered
    bool heightmapFlag_;            ///< flag whether a height map will be rendered instead of a surface

    // scatter plot related settings:
    float minGlyphSize_;            ///< minimum glyph size
    float maxGlyphSize_;            ///< maximum glyph size
    GlyphStyle glyphStyle_;         ///< glyph style
    std::string texturePath_;       ///< texture path
    bool useTextureFlag_;           ///< flag whether to use textured glyphs

};

} // namespace voreen

#endif // VRN_PLOTENTITYSETTINGS_H
