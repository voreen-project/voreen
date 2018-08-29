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

#include "plotentitysettings.h"

#include "voreen/core/io/serialization/serialization.h"
#include "tgt/filesystem.h"

#include <vector>

namespace voreen {

PlotEntitySettings::PlotEntitySettings():
    entity_(NONE)
    , mainColumnIndex_(-1)
    , candleTopColumnIndex_(-1)
    , candleBottomColumnIndex_(-1)
    , stickTopColumnIndex_(-1)
    , stickBottomColumnIndex_(-1)
    , optionalColumnIndex_(-1)
    , secondOptionalColumnIndex_(-1)
    , firstColor_(tgt::Color(0.f))
    , secondColor_(tgt::Color(0.f))
    , colorMap_(ColorMap::createColorMap(0))
    , lineStyle_(CONTINUOUS)
    , splineFlag_(false)
    , errorbarFlag_(false)
    , candleStickFlag_(false)
    , wireOnlyFlag_(false)
    , heightmapFlag_(false)
    , minGlyphSize_(1)
    , maxGlyphSize_(1)
    , glyphStyle_(POINT)
    , texturePath_("")
    , useTextureFlag_(false)
{}

PlotEntitySettings::PlotEntitySettings(int lineCI, tgt::Color lineColor, LineStyle style,
        bool splineFlag, int errorCI, tgt::Color errorColor,  bool errorbarFlag):
    entity_(LINE)
    , mainColumnIndex_(lineCI)
    , candleTopColumnIndex_(lineCI)
    , candleBottomColumnIndex_(lineCI)
    , stickTopColumnIndex_(lineCI)
    , stickBottomColumnIndex_(lineCI)
    , optionalColumnIndex_(errorCI)
    , secondOptionalColumnIndex_(-1)
    , firstColor_(lineColor)
    , secondColor_(errorColor)
    , colorMap_(ColorMap::createColorMap(0))
    , lineStyle_(style)
    , splineFlag_(splineFlag)
    , errorbarFlag_(errorbarFlag)
    , candleStickFlag_(false)
    , wireOnlyFlag_(false)
    , heightmapFlag_(false)
    , minGlyphSize_(1)
    , maxGlyphSize_(1)
    , glyphStyle_(POINT)
    , texturePath_("")
    , useTextureFlag_(false)
    {}

PlotEntitySettings::PlotEntitySettings(int barsCI, tgt::Color barColor):
    entity_(BAR)
    , mainColumnIndex_(barsCI)
    , candleTopColumnIndex_(0)
    , candleBottomColumnIndex_(0)
    , stickTopColumnIndex_(0)
    , stickBottomColumnIndex_(0)
    , optionalColumnIndex_(-1)
    , secondOptionalColumnIndex_(-1)
    , firstColor_(barColor)
    , secondColor_(barColor)
    , colorMap_(ColorMap::createColorMap(0))
    , lineStyle_(CONTINUOUS)
    , splineFlag_(false)
    , errorbarFlag_(false)
    , candleStickFlag_(false)
    , wireOnlyFlag_(false)
    , heightmapFlag_(false)
    , minGlyphSize_(1)
    , maxGlyphSize_(1)
    , glyphStyle_(POINT)
    , texturePath_("")
    , useTextureFlag_(false)
    {}


PlotEntitySettings::PlotEntitySettings(int surfaceCI, tgt::Color surfaceColor, bool wireFlag, tgt::Color wireColor,
    ColorMap cm, bool heightmapFlag, int optionalCI):
    entity_(SURFACE)
    , mainColumnIndex_(surfaceCI)
    , candleTopColumnIndex_(0)
    , candleBottomColumnIndex_(0)
    , stickTopColumnIndex_(0)
    , stickBottomColumnIndex_(0)
    , optionalColumnIndex_(optionalCI)
    , secondOptionalColumnIndex_(-1)
    , firstColor_(surfaceColor)
    , secondColor_(wireColor)
    , colorMap_(cm)
    , lineStyle_(CONTINUOUS)
    , splineFlag_(false)
    , errorbarFlag_(false)
    , candleStickFlag_(false)
    , wireOnlyFlag_(wireFlag)
    , heightmapFlag_(heightmapFlag)
    , minGlyphSize_(1)
    , maxGlyphSize_(1)
    , glyphStyle_(POINT)
    , texturePath_("")
    , useTextureFlag_(false)
    {}


PlotEntitySettings::PlotEntitySettings(int scatterCI, tgt::Color color, float glyphSizeMin ,float glyphSizeMax, GlyphStyle style):
    entity_(SCATTER)
    , mainColumnIndex_(scatterCI)
    , candleTopColumnIndex_(0)
    , candleBottomColumnIndex_(0)
    , stickTopColumnIndex_(0)
    , stickBottomColumnIndex_(0)
    , optionalColumnIndex_(-1)
    , secondOptionalColumnIndex_(-1)
    , firstColor_(color)
    , secondColor_(color)
    , colorMap_(ColorMap::createColorMap(0))
    , lineStyle_(CONTINUOUS)
    , splineFlag_(false)
    , errorbarFlag_(false)
    , candleStickFlag_(false)
    , wireOnlyFlag_(false)
    , heightmapFlag_(false)
    , minGlyphSize_(glyphSizeMin)
    , maxGlyphSize_(glyphSizeMax)
    , glyphStyle_(style)
    , texturePath_("")
    , useTextureFlag_(false)
    {}

void PlotEntitySettings::serialize(Serializer& s) const {
    s.serialize("entity", static_cast<int>(entity_));
    s.serialize("colorMap", colorMap_);
    s.serialize("mainColumnIndex", mainColumnIndex_);
    s.serialize("candleTopColumnIndex", candleTopColumnIndex_);
    s.serialize("candleBottomColumnIndex", candleBottomColumnIndex_);
    s.serialize("stickTopColumnIndex", stickTopColumnIndex_);
    s.serialize("stickBottomColumnIndex", stickBottomColumnIndex_);
    s.serialize("optionalCI", optionalColumnIndex_);
    s.serialize("secondOptionalCI", secondOptionalColumnIndex_);
    s.serialize("firstColor", firstColor_);
    s.serialize("secondColor", secondColor_);
    s.serialize("lineStyle", static_cast<int>(lineStyle_));
    s.serialize("splineFlag", splineFlag_);
    s.serialize("errorbarFlag", errorbarFlag_);
    s.serialize("wireOnlyFlag", wireOnlyFlag_);
    s.serialize("heightmapFlag", heightmapFlag_);
    s.serialize("candleStickFlag", candleStickFlag_);
    s.serialize("minGlyphSize", minGlyphSize_);
    s.serialize("maxGlyphSize", maxGlyphSize_);
    s.serialize("glyphStyle", static_cast<int>(glyphStyle_));

    //see FileDialogProperty
    // convert path to an relative one with respect to the document's path
    std::string path = texturePath_;
    if (!path.empty() && !s.getDocumentPath().empty())
        path = tgt::FileSystem::relativePath(path, tgt::FileSystem::dirName(s.getDocumentPath()));

    // cleanup path: replace backslashes
    std::string::size_type pos = path.find("\\");
    while (pos != std::string::npos) {
        path[pos] = '/';
        pos = path.find("\\");
    }
    s.serialize("texturePath", path);
    s.serialize("useTextureFlag", useTextureFlag_);
}

void PlotEntitySettings::deserialize(Deserializer& s) {
    int value;
    s.deserialize("entity", value);
    entity_ = static_cast<Entity>(value);
    s.deserialize("colorMap", colorMap_);
    s.deserialize("mainColumnIndex", mainColumnIndex_);
    s.deserialize("candleTopColumnIndex", candleTopColumnIndex_);
    s.deserialize("candleBottomColumnIndex", candleBottomColumnIndex_);
    s.deserialize("stickTopColumnIndex", stickTopColumnIndex_);
    s.deserialize("stickBottomColumnIndex", stickBottomColumnIndex_);
    s.deserialize("optionalCI", optionalColumnIndex_);
    s.deserialize("secondOptionalCI", secondOptionalColumnIndex_);
    s.deserialize("firstColor", firstColor_);
    s.deserialize("secondColor", secondColor_);
    s.deserialize("lineStyle", value);
    lineStyle_ = static_cast<LineStyle>(value);
    s.deserialize("splineFlag", splineFlag_);
    s.deserialize("errorbarFlag", errorbarFlag_);
    s.deserialize("wireOnlyFlag", wireOnlyFlag_);
    s.deserialize("heightmapFlag", heightmapFlag_);
    s.deserialize("candleStickFlag", candleStickFlag_);
    s.deserialize("minGlyphSize", minGlyphSize_);
    s.deserialize("maxGlyphSize", maxGlyphSize_);
    s.deserialize("glyphStyle", value);
    glyphStyle_ = static_cast<GlyphStyle>(value);
    //see FileDialogProperty
    s.deserialize("texturePath", texturePath_);
    // convert path relative to the document's path to an absolute one
    if (!texturePath_.empty() && !s.getDocumentPath().empty())
        texturePath_ = tgt::FileSystem::absolutePath(
            tgt::FileSystem::dirName(s.getDocumentPath()) + "/" + texturePath_);
    s.deserialize("useTextureFlag", useTextureFlag_);
}

bool PlotEntitySettings::fitsTo(const PlotData* data) const {
    int cc = data->getColumnCount();
    if (entity_ == LINE) {
        if (mainColumnIndex_ < cc && candleTopColumnIndex_ < cc
            && candleBottomColumnIndex_ < cc && stickTopColumnIndex_ < cc
            && stickBottomColumnIndex_ < cc && optionalColumnIndex_ < cc
            && data->getColumnType(mainColumnIndex_)== PlotBase::NUMBER
            && (optionalColumnIndex_ == -1 || data->getColumnType(optionalColumnIndex_)== PlotBase::NUMBER)
            && data->getColumnType(candleTopColumnIndex_)== PlotBase::NUMBER
            && data->getColumnType(candleBottomColumnIndex_)== PlotBase::NUMBER
            && data->getColumnType(stickTopColumnIndex_)== PlotBase::NUMBER
            && data->getColumnType(stickBottomColumnIndex_)== PlotBase::NUMBER)

            return true;
        else
            return false;
    }
    else if (entity_ == BAR) {
        if (mainColumnIndex_ < cc && data->getColumnType(mainColumnIndex_)== PlotBase::NUMBER)
            return true;
        else
            return false;
    }
    else if (entity_ == SURFACE) {
        if (mainColumnIndex_ < cc && optionalColumnIndex_ < cc
            && data->getColumnType(mainColumnIndex_)== PlotBase::NUMBER
            && (optionalColumnIndex_ == -1 || data->getColumnType(optionalColumnIndex_)== PlotBase::NUMBER))
            return true;
        else
            return false;
    }
    else { // SCATTER
        if (mainColumnIndex_ < cc && candleTopColumnIndex_ < cc
            && data->getColumnType(mainColumnIndex_)== PlotBase::NUMBER
            && (optionalColumnIndex_ == -1 || data->getColumnType(optionalColumnIndex_)== PlotBase::NUMBER)
            && (secondOptionalColumnIndex_ == -1 || data->getColumnType(secondOptionalColumnIndex_)== PlotBase::NUMBER))
            return true;
        else
            return false;
    }
}

bool PlotEntitySettings::getCandleStickFlag() const {
    return candleStickFlag_;
}

void PlotEntitySettings::setCandleStickFlag(bool value) {
    candleStickFlag_ = value;
}

PlotEntitySettings::Entity PlotEntitySettings::getEntity() const {
    return entity_;
}

void PlotEntitySettings::setEntity(Entity entity) {
    entity_ = entity;
}

int PlotEntitySettings::getMainColumnIndex() const {
    return mainColumnIndex_;
}

void PlotEntitySettings::setMainColumnIndex(int value) {
    mainColumnIndex_ = value;
}

int PlotEntitySettings::getCandleTopColumnIndex() const {
    return candleTopColumnIndex_;
}
void PlotEntitySettings::setCandleTopColumnIndex(int value) {
    candleTopColumnIndex_ = value;
}

int PlotEntitySettings::getCandleBottomColumnIndex() const {
    return candleBottomColumnIndex_;
}
void PlotEntitySettings::setCandleBottomColumnIndex(int value) {
    candleBottomColumnIndex_ = value;
}

int PlotEntitySettings::getStickTopColumnIndex() const {
    return stickTopColumnIndex_;
}
void PlotEntitySettings::setStickTopColumnIndex(int value) {
    stickTopColumnIndex_ = value;
}

int PlotEntitySettings::getStickBottomColumnIndex() const {
    return stickBottomColumnIndex_;
}
void PlotEntitySettings::setStickBottomColumnIndex(int value) {
    stickBottomColumnIndex_ = value;
}

int PlotEntitySettings::getOptionalColumnIndex() const {
    return optionalColumnIndex_;

}
void PlotEntitySettings::setOptionalColumnIndex(int value) {
    optionalColumnIndex_ = value;
}

tgt::Color PlotEntitySettings::getFirstColor() const {
    return firstColor_;
}
void PlotEntitySettings::setFirstColor(tgt::Color value) {
    firstColor_ = value;
}

tgt::Color PlotEntitySettings::getSecondColor() const {
    return secondColor_;
}
void PlotEntitySettings::setSecondColor(tgt::Color value) {
    secondColor_ = value;
}

PlotEntitySettings::LineStyle PlotEntitySettings::getLineStyle() const {
    return lineStyle_;
}
void PlotEntitySettings::setLineStyle(PlotEntitySettings::LineStyle value) {
    lineStyle_ = value;
}

bool PlotEntitySettings::getSplineFlag() const {
    return splineFlag_;
}

void PlotEntitySettings::setSplineFlag(bool value) {
    splineFlag_ = value;
}

bool PlotEntitySettings::getErrorbarFlag() const {
    return errorbarFlag_;
}
void PlotEntitySettings::setErrorbarFlag(bool value) {
    errorbarFlag_ = value;
}

bool PlotEntitySettings::getWireOnlyFlag() const {
    return wireOnlyFlag_;
}

void PlotEntitySettings::setWireOnlyFlag(bool value) {
    wireOnlyFlag_ = value;
}

bool PlotEntitySettings::getHeightmapFlag() const {
    return heightmapFlag_;
}

void PlotEntitySettings::setHeightmapFlag(bool value) {
    heightmapFlag_ = value;
}

ColorMap PlotEntitySettings::getColorMap() const {
    return colorMap_;
}

void PlotEntitySettings::setColorMap(ColorMap colorMap) {
    colorMap_ = colorMap;
}

float PlotEntitySettings::getMinGlyphSize() const {
    return minGlyphSize_;
}

void PlotEntitySettings::setMinGlyphSize(float value) {
    minGlyphSize_ = value;
}

float PlotEntitySettings::getMaxGlyphSize() const {
    return maxGlyphSize_;
}

void PlotEntitySettings::setMaxGlyphSize(float value) {
    maxGlyphSize_ = value;
}

PlotEntitySettings::GlyphStyle PlotEntitySettings::getGlyphStyle() const {
    return glyphStyle_;
}

void PlotEntitySettings::setGlyphStyle(PlotEntitySettings::GlyphStyle value) {
    glyphStyle_ = value;
}

int PlotEntitySettings::getSecondOptionalColumnIndex() const {
    return secondOptionalColumnIndex_;
}

void PlotEntitySettings::setSecondOptionalColumnIndex( int value ) {
    secondOptionalColumnIndex_ = value;
}

std::string PlotEntitySettings::getTexturePath() const {
    return texturePath_;
}

void PlotEntitySettings::setTexturePath(std::string path) {
    texturePath_ = path;
}

bool PlotEntitySettings::getUseTextureFlag() const {
    return useTextureFlag_;
}

void PlotEntitySettings::setUseTextureFlag(bool value) {
    useTextureFlag_ = value;
}

bool PlotEntitySettings::operator==(const PlotEntitySettings& rhs) const {
    if (entity_ != rhs.getEntity())
        return false;

    if (entity_ == NONE)
        return true;
    else if (entity_ == LINE) {
        if (mainColumnIndex_ == rhs.getMainColumnIndex() &&
            candleTopColumnIndex_ == rhs.getCandleTopColumnIndex() &&
            candleBottomColumnIndex_ == rhs.getCandleBottomColumnIndex() &&
            stickTopColumnIndex_ == rhs.getStickTopColumnIndex() &&
            stickBottomColumnIndex_ == rhs.getStickBottomColumnIndex() &&
            firstColor_ == rhs.getFirstColor() &&
            secondColor_ == rhs.getSecondColor() &&
            lineStyle_ == rhs.getLineStyle() &&
            splineFlag_ == rhs.getSplineFlag() &&
            errorbarFlag_ == rhs.getErrorbarFlag() &&
            candleStickFlag_ == rhs.getCandleStickFlag())
            return true;
    }
    else if (entity_ == BAR) {
        if (mainColumnIndex_ == rhs.getMainColumnIndex() &&
            firstColor_ == rhs.getFirstColor())
            return true;
    }
    else if (entity_ == SURFACE) {
        if (mainColumnIndex_ == rhs.getMainColumnIndex() &&
            optionalColumnIndex_ == rhs.getOptionalColumnIndex() &&
            firstColor_ == rhs.getFirstColor() &&
            secondColor_ == rhs.getSecondColor() &&
            wireOnlyFlag_ == rhs.getWireOnlyFlag() &&
            heightmapFlag_ == rhs.getHeightmapFlag() &&
            colorMap_ == rhs.getColorMap())
            return true;
    }
    else if (entity_ == SCATTER) {
        if (mainColumnIndex_ == rhs.getMainColumnIndex() &&
            optionalColumnIndex_ == rhs.getOptionalColumnIndex() &&
            secondOptionalColumnIndex_ == rhs.getSecondOptionalColumnIndex() &&
            firstColor_ == rhs.getFirstColor() &&
            colorMap_ == rhs.getColorMap() &&
            minGlyphSize_ == rhs.getMinGlyphSize() &&
            maxGlyphSize_ == rhs.getMaxGlyphSize() &&
            glyphStyle_ == rhs.getGlyphStyle() &&
            texturePath_ == rhs.getTexturePath() &&
            useTextureFlag_ == rhs.getUseTextureFlag())
            return true;
    }
    return false;

}

} // namespace voreen
