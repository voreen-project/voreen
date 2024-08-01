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

#include "plotentitiesproperty.h"

namespace voreen {

PlotEntitiesProperty::PlotEntitiesProperty(const std::string& id, const std::string& guiText, PlotEntitySettings::Entity entities,
       std::vector<PlotEntitySettings> value, Processor::InvalidationLevel invalidationLevel)
       : TemplateProperty<std::vector<PlotEntitySettings> >(id, guiText, value, invalidationLevel)
    , entities_(entities)
    , xColumnIndex_(0)
    , yColumnIndex_(-1)
    , colorMap_(ColorMap::createColorMap(0))
    , dataValidFlag_(false)
    , loadStrategy_(LS_ALL)
{}

PlotEntitiesProperty::PlotEntitiesProperty()
    : colorMap_(ColorMap::createColorMap(0))
{}

Property* PlotEntitiesProperty::create() const {
    return new PlotEntitiesProperty();
}

void PlotEntitiesProperty::serialize(Serializer& s) const {
    Property::serialize(s);

    s.serialize("entities", static_cast<int>(entities_));
    s.serialize("xCI", xColumnIndex_);
    s.serialize("yCI", yColumnIndex_);
    s.serialize("colorMap", colorMap_);
    s.serialize("plotEntitySettingsVector", value_, "plotEntitySettings");
    s.serialize("loadStrategy",static_cast<int>(loadStrategy_));
}

void PlotEntitiesProperty::deserialize(Deserializer& s) {
    Property::deserialize(s);

    int value, loadStrategy;
    s.deserialize("entities", value);
    entities_ = static_cast<PlotEntitySettings::Entity>(value);
    s.deserialize("xCI", xColumnIndex_);
    s.deserialize("yCI", yColumnIndex_);
    s.deserialize("colorMap", colorMap_);
    s.deserialize("plotEntitySettingsVector", value_, "plotEntitySettings");
    s.deserialize("loadStrategy",loadStrategy);
    loadStrategy_ = static_cast<PlotEntitiesProperty::loadStrategy>(loadStrategy);
}

PlotEntitySettings::Entity PlotEntitiesProperty::getEntities() const {
    return entities_;
}

void PlotEntitiesProperty::setEntities(PlotEntitySettings::Entity entities) {
    entities_ = entities;
}

int PlotEntitiesProperty::getXColumnIndex() const {
    return xColumnIndex_;
}

bool PlotEntitiesProperty::setXColumnIndex(int index) {
    if (dataEmpty() || index < 0 || index >= data_->getColumnCount())
        return false;
    else {
        xColumnIndex_ = index;
        invalidate();
        return true;
    }
}

PlotEntitiesProperty::loadStrategy PlotEntitiesProperty::getLoadStrategy() const {
    return loadStrategy_;
}

void PlotEntitiesProperty::setLoadStrategy(loadStrategy ls) {
    loadStrategy_ = ls;
}

int PlotEntitiesProperty::getYColumnIndex() const {
    if (entities_ == PlotEntitySettings::SCATTER || entities_ == PlotEntitySettings::SURFACE)
        return yColumnIndex_;
    else
        throw VoreenException("y value exist only for SCATTER and SURFACE");
}

bool PlotEntitiesProperty::setYColumnIndex(int index) {
    if (dataEmpty() || index >= data_->getColumnCount())
        return false;
    else if (entities_ != PlotEntitySettings::SURFACE && entities_ != PlotEntitySettings::SCATTER)
        return false;
    else {
        yColumnIndex_ = index;
        invalidate();
        return true;
    }
}

std::vector<int> PlotEntitiesProperty::getDataColumnIndices() const {
    std::vector<int> toReturn;
    for (std::vector<PlotEntitySettings>::const_iterator it = value_.begin(); it != value_.end(); ++it) {
        if (it->getCandleStickFlag()) {
            toReturn.push_back(it->getCandleBottomColumnIndex());
            toReturn.push_back(it->getCandleTopColumnIndex());
            toReturn.push_back(it->getStickBottomColumnIndex());
            toReturn.push_back(it->getStickTopColumnIndex());
        }
        else
            toReturn.push_back(it->getMainColumnIndex());
    }
    return toReturn;
}

const ColorMap& PlotEntitiesProperty::getColorMap() const {
    return colorMap_;
}

void PlotEntitiesProperty::setColorMap(ColorMap cm) {
    colorMap_ = cm;
    invalidate();
}

const PlotData* PlotEntitiesProperty::getPlotData() const {
    return data_;
}

void PlotEntitiesProperty::setPlotData(const PlotData* data) {
    data_ = data;
    dataValidFlag_ = false;

    int possibleXIndex = 0;
    int possibleYIndex = -1;

    if (!dataEmpty()) {
        switch (entities_) {
            case PlotEntitySettings::LINE:
            case PlotEntitySettings::BAR: {
                //check if there is a possible x axis
                bool xAxisFound = false;
                //check if there is another number column as y axis
                bool dataAxisFound = false;
                for (int i = 0; i < data_->getColumnCount(); ++i) {
                    if (possibleXAxis(i) && !xAxisFound) {
                        xAxisFound = true;
                        possibleXIndex = i;
                        continue;
                    }
                    if (data_->getColumnType(i) == PlotBase::NUMBER)
                        dataAxisFound = true;
                }
                dataValidFlag_ = (xAxisFound && dataAxisFound);
                break;
            }
            case PlotEntitySettings::SURFACE: {
                bool xAxisFound = false;
                bool yAxisFound = false;
                bool dataAxisFound = false;
                for (int i = 0; i < data_->getColumnCount(); ++i) {
                    if (possibleXAxis(i) && !xAxisFound) {
                        xAxisFound = true;
                        possibleXIndex = i;
                        continue;
                    }
                    if (possibleYAxis(i) && !yAxisFound) {
                        yAxisFound = true;
                        possibleYIndex = i;
                        continue;
                    }
                    if (data_->getColumnType(i) == PlotBase::NUMBER)
                        dataAxisFound = true;
                }
                dataValidFlag_ = (xAxisFound && yAxisFound && dataAxisFound);
                break;
            }
            case PlotEntitySettings::SCATTER: {
                bool xAxisFound = false;
                bool dataAxisFound = false;
                for (int i = 0; i < data_->getColumnCount(); ++i) {
                    if (possibleXAxis(i) && !xAxisFound) {
                        xAxisFound = true;
                        possibleXIndex = i;
                        continue;
                    }
                    if (data_->getColumnType(i) == PlotBase::NUMBER)
                        dataAxisFound = true;
                }
                dataValidFlag_ = (xAxisFound && dataAxisFound);
                break;
            }
            default :
                break;
        }
    }

    //possible reset
    if (dataEmpty() || !possibleXAxis(xColumnIndex_) || !possibleYAxis(yColumnIndex_)) {
        xColumnIndex_ = possibleXIndex;
        yColumnIndex_ = possibleYIndex;
        value_.clear();
        if (possibleXAxis(xColumnIndex_) && possibleYAxis(yColumnIndex_)){
            value_ = createAllEntitySettings();
        }
    }
    else {
        switch(loadStrategy_)
        {
        case LS_ALL:
            xColumnIndex_ = possibleXIndex;
            yColumnIndex_ = possibleYIndex;
            value_.clear();
            if (possibleXAxis(xColumnIndex_) && possibleYAxis(yColumnIndex_)){
                value_ = createAllEntitySettings();
            }
            break;
        case LS_NON:
        case LS_NEW:
            bool reset = false;
            for (size_t i = 0; i < value_.size(); ++i) {
                if (!value_.at(i).fitsTo(data_)) {
                    xColumnIndex_ = possibleXIndex;
                    yColumnIndex_ = possibleYIndex;
                    value_.clear();
                    if (possibleXAxis(xColumnIndex_) && possibleYAxis(yColumnIndex_)){
                        value_ = createAllEntitySettings();
                    }
                    reset = true;
                    break;
                }
            }

            // adjust colors to color hints:
            if(!reset) {
                for (size_t i=0; i<value_.size(); ++i) {
                    int ci = value_.at(i).getMainColumnIndex();
                    if (data_->hasColumnColorHint(ci)) {
                        value_[i].setFirstColor(data_->getColumnColorHint(ci));
                    }
                }
            }

            if(loadStrategy_ == LS_NON)
                break;

            //FIXME: if size is not well
            if(value_.size() < static_cast<size_t>(data_->getDataColumnCount())){
                value_ = addEntitySettings(); //< function also not well
            }
        }

    }
    invalidate();
}

void PlotEntitiesProperty::setPlotEntitySettings(PlotEntitySettings settings, int index) {
    if (entities_ == settings.getEntity()) {
        if (!dataEmpty()) {
            //new entity
            if (static_cast<size_t>(index) >= value_.size())
                value_.push_back(settings);
            else {
                value_.erase(value_.begin()+index);
                value_.insert(value_.begin()+index, settings);
            }
            invalidate();
        }
    }
    else
        throw VoreenException("PlotEntitiesProperty::setPlotEntity: entities do not fit");
}

void PlotEntitiesProperty::deletePlotEntitySettings(int index) {
    if (!dataEmpty()) {
        if (static_cast<size_t>(index) < value_.size())
            value_.erase(value_.begin()+index);
        invalidate();
    }
}

bool PlotEntitiesProperty::dataEmpty() const {
    if (data_ && data_->getColumnCount() > 0 && !data_->rowsEmpty())
        return false;
    else
        return true;
}

bool PlotEntitiesProperty::dataValid() const {
    return dataValidFlag_;
}

void PlotEntitiesProperty::applyColormap() {
    if (!dataEmpty()) {
        ColorMap::GeneratingIterator cmit = colorMap_.getGeneratingIterator();
        for (size_t i = 0; i < value_.size(); ++i) {
            PlotEntitySettings es(value_.at(i));
            if(!(data_ && ((int)i < data_->getColumnCount()) && !data_->hasColumnColorHint((int)i))) {
                es.setFirstColor(*cmit);
                ++cmit;
                if (es.getEntity() == PlotEntitySettings::LINE) {
                    es.setSecondColor(*cmit);
                    //++cmit;
                }
                else if (es.getEntity() == PlotEntitySettings::SURFACE) {
                    es.setSecondColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
                }
            }
            setPlotEntitySettings(es, static_cast<int>(i));
        }
    }
}

PlotEntitySettings PlotEntitiesProperty::createEntitySettings() const {
    if (!dataValidFlag_) {
        return PlotEntitySettings();
    }
    else {
        PlotEntitySettings es;
        ColorMap::GeneratingIterator it = ColorMap::GeneratingIterator(&colorMap_, static_cast<int>(value_.size()));
        //find number column index
        int ci;
        for (ci = 0; ci < data_->getColumnCount(); ++ci) {
            if (data_->getColumnType(ci)==PlotBase::NUMBER
                && xColumnIndex_ != ci
                && (entities_ == PlotEntitySettings::LINE || entities_ == PlotEntitySettings::BAR || yColumnIndex_ != ci))
                break;
        }
        tgt::Color c = *it;
        if(data_->hasColumnColorHint(ci))
            c = data_->getColumnColorHint(ci);
        switch (entities_) {
            case PlotEntitySettings::LINE:
                es = PlotEntitySettings(ci,c,PlotEntitySettings::CONTINUOUS, false, -1,c,false);
                break;
            case PlotEntitySettings::BAR:
                es = PlotEntitySettings(ci,c);
                break;
            case PlotEntitySettings::SURFACE:
                es = PlotEntitySettings(ci,c, false, tgt::Color(0.f, 0.f, 0.f, 1.f), ColorMap::createColorMap(0), false, -1);
                break;
            case PlotEntitySettings::SCATTER:
                es = PlotEntitySettings(ci, c, 1, 2, PlotEntitySettings::POINT);
                break;
            default :
                es = PlotEntitySettings();
                throw VoreenException("PlotEntitiesProperty::createEntitySettings: this shouldn't be reached");
                break;
        }
        return es;
    }
}

std::vector<PlotEntitySettings> PlotEntitiesProperty::createAllEntitySettings() const {
    if (!dataValidFlag_) {
        return std::vector<PlotEntitySettings>();
    }
    else {
        std::vector<PlotEntitySettings> all;
        ColorMap::GeneratingIterator it = ColorMap::GeneratingIterator(&colorMap_,
            static_cast<int>(value_.size()));
        for (int ci = 0; ci < data_->getColumnCount(); ++ci) {
            PlotEntitySettings es;
            //find number column index
            if (data_->getColumnType(ci)==PlotBase::NUMBER
                    && xColumnIndex_ != ci
                    && (entities_ == PlotEntitySettings::LINE || entities_ == PlotEntitySettings::BAR || yColumnIndex_ != ci)) {
                tgt::Color c = *it;
                if(data_->hasColumnColorHint(ci))
                    c = data_->getColumnColorHint(ci);

                switch (entities_) {
                    case PlotEntitySettings::LINE:
                        es = PlotEntitySettings(ci,c,PlotEntitySettings::CONTINUOUS, false, -1,c,false);
                        break;
                    case PlotEntitySettings::BAR:
                        es = PlotEntitySettings(ci,c);
                        break;
                    case PlotEntitySettings::SURFACE:
                        es = PlotEntitySettings(ci,c, false, tgt::Color(0.f, 0.f, 0.f, 1.f), ColorMap::createColorMap(0), false, -1);
                        break;
                    case PlotEntitySettings::SCATTER:
                        es = PlotEntitySettings(ci, c, 1, 2, PlotEntitySettings::POINT);
                        break;
                    default :
                        es = PlotEntitySettings();
                        throw VoreenException("PlotEntitiesProperty::createEntitySettings: this shouldn't be reached");
                        break;
                }
                all.push_back(es);
                if(!data_->hasColumnColorHint(ci))
                    ++it;
            }
        }
        return all;
    }
}

std::vector<PlotEntitySettings> PlotEntitiesProperty::addEntitySettings() const {
    if (!dataValidFlag_) {
        return std::vector<PlotEntitySettings>();
    }
    else {
        std::vector<PlotEntitySettings> all(value_);
        ColorMap::GeneratingIterator it = ColorMap::GeneratingIterator(&colorMap_,
            static_cast<int>(value_.size()));
        for (int ci = static_cast<int>(value_.size()); ci < data_->getColumnCount(); ++ci) {
            PlotEntitySettings es;
            //find number column index
            if (data_->getColumnType(ci)==PlotBase::NUMBER
                    && xColumnIndex_ != ci
                    && (entities_ == PlotEntitySettings::LINE || entities_ == PlotEntitySettings::BAR || yColumnIndex_ != ci)) {
                switch (entities_) {
                    case PlotEntitySettings::LINE:
                        es = PlotEntitySettings(ci,*it,PlotEntitySettings::CONTINUOUS, false, -1,*it,false);
                        break;
                    case PlotEntitySettings::BAR:
                        es = PlotEntitySettings(ci,*it);
                        break;
                    case PlotEntitySettings::SURFACE:
                        es = PlotEntitySettings(ci,*it, false, tgt::Color(0.f, 0.f, 0.f, 1.f), ColorMap::createColorMap(0), false, -1);
                        break;
                    case PlotEntitySettings::SCATTER:
                        es = PlotEntitySettings(ci, *it, 1, 2, PlotEntitySettings::POINT);
                        break;
                    default :
                        es = PlotEntitySettings();
                        throw VoreenException("PlotEntitiesProperty::createEntitySettings: this shouldn't be reached");
                        break;
                }
                all.push_back(es);
                ++it;
            }
        }
        return all;
    }
}

bool PlotEntitiesProperty::possibleXAxis(int index) {
    //check if data empty or index out of bounds
    if (dataEmpty() || index < 0 || index >= data_->getColumnCount())
        return false;

    //for LINE and BAR, label x axis is possible
    if (data_->getColumnType(index) == PlotBase::STRING) {
        if (entities_ == PlotEntitySettings::LINE || entities_ == PlotEntitySettings::BAR)
            return true;
        else
            return false;
    }
    //number column is good enough (but for LINE it must be a key column)
    if (data_->getColumnType(index) == PlotBase::NUMBER && (entities_ != PlotEntitySettings::LINE || index <= data_->getKeyColumnCount()))
        return true;
    else
        return false;
}

bool PlotEntitiesProperty::possibleYAxis(int index) {
    //check if data empty or index out of bounds
    if (dataEmpty() || index < -1 || index >= data_->getColumnCount())
        return false;

    //-1 isn't possible for SURFACE otherwise is it possible
    if (index == -1) {
        if (entities_ == PlotEntitySettings::SURFACE)
            return false;
        else
            return true;
    }

    if ((entities_ == PlotEntitySettings::SURFACE || entities_ == PlotEntitySettings::SCATTER)
        && data_->getColumnType(index) == PlotBase::NUMBER)
        return true;
    else
        return false;
}

/*
up to now, we do not support linking of this property. If it will be implemented, this method must check
whether the structure of the PlotData match. This is the reason why we need to overwrite this
method. To implement the linking, one has to create a new class storing vector<PlotEntitySettings> and
x axis column index etc.
void PlotEntitiesProperty::executeLinks(const std::vector<PlotEntitySettings>& prevValue, const std::vector<PlotEntitySettings>& curValue) {

    if (links_.empty())
        return;

    // check if all properties work with the same PlotData
    for (std::vector<PropertyLink*>::iterator it = links_.begin(); it != links_.end(); ++it) {
        PlotEntitiesProperty* pep = dynamic_cast<PlotEntitiesProperty*>((*it)->getDestinationProperty());
        if (!pep || !dataValidFlag_ || !pep->dataValid() || !data_->compareStructure(pep->getPlotData()))
            return;
    }

    // pass change data object to links
    for (std::vector<PropertyLink*>::iterator it = links_.begin(); it != links_.end(); it++) {
        try {
            (*it)->onChange();
        }
        catch (const VoreenException& e) {
            LERRORC("voreen.TemplateProperty", "executeLinks(): " << e.what());
        }
    }
}*/

} // namespace voreen
