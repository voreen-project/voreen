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

#ifndef VRN_PLOTENTITIESPROPERTY_H
#define VRN_PLOTENTITIESPROPERTY_H

#include "voreen/core/properties/templateproperty.h"
#include "../datastructures/colormap.h"
#include "../datastructures/plotdata.h"
#include "../datastructures/plotentitysettings.h"

#include <vector>

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API TemplateProperty<std::vector<PlotEntitySettings> >;
#endif
/**
 * Property encapsulating a vector of PlotEntitySettings.
 */
class VRN_CORE_API PlotEntitiesProperty : public TemplateProperty<std::vector<PlotEntitySettings> > {
public:
     ///load strategy
    enum loadStrategy {
        LS_NON,             ///< load no data
        LS_NEW,             ///< load only new data column
        LS_ALL              ///< load all data column
    };

    PlotEntitiesProperty(const std::string& id, const std::string& guiText, PlotEntitySettings::Entity entities,
       std::vector<PlotEntitySettings> value = std::vector<PlotEntitySettings>(), Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT);
    PlotEntitiesProperty();

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "PlotEntitiesProperty"; }
    virtual std::string getTypeDescription() const { return "PlotEntities"; }

    /**
     * Executes all member actions that belong to the property. Generally the owner of the
     * property is invalidated.
     */
    //void notifyChange();

    /**
     * @see Property::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Property::deserialize
     */
    virtual void deserialize(Deserializer& s);


    PlotEntitySettings::Entity getEntities() const;
    void setEntities(PlotEntitySettings::Entity entities);

    ///current selected x axis
    int getXColumnIndex() const;
    bool setXColumnIndex(int index);

    //call only if entities == SURFACE or SCATTER
    ///current selected y axis
    int getYColumnIndex() const;
    bool setYColumnIndex(int index);

    /**
     * Returns a vector of the main indices of all held PlotEntities
     **/
    std::vector<int> getDataColumnIndices() const;

    const ColorMap& getColorMap() const;
    void setColorMap(ColorMap cm);

    const PlotData* getPlotData() const;
    void setPlotData(const PlotData* data);

    void setPlotEntitySettings(PlotEntitySettings settings, int index);
    void deletePlotEntitySettings(int index);

    loadStrategy getLoadStrategy() const;
    void setLoadStrategy(loadStrategy ls);

    bool dataEmpty() const;

    /// returns the dataValidFlag_
    bool dataValid() const;

    ///applies the colormap to the entities
    void applyColormap();

    ///creates entity settings
    PlotEntitySettings createEntitySettings() const;

    ///creates all entity settings
    std::vector<PlotEntitySettings> createAllEntitySettings() const;

    ///add entity settings
    std::vector<PlotEntitySettings> addEntitySettings() const;

private:
    /// executes links if all properties use the same PlotData
    //virtual void executeLinks(const std::vector<PlotEntitySettings>& prevValue, const std::vector<PlotEntitySettings>& curValue);

    ///returns true, if the column of data_ with index index can represent a x axis
    bool possibleXAxis(int index);
    ///returns true, if the column of data_ with index index can represent a y axis
    bool possibleYAxis(int index);

    /// type of the entities
    PlotEntitySettings::Entity entities_;

    int xColumnIndex_;                                ///< column index of x axis (ALL)
    int yColumnIndex_;                                ///< column index of y axis (SURFACE and SCATTER)
    ColorMap colorMap_;

    bool dataValidFlag_;     ///< false if there is not enough data to plot something

    const PlotData* data_;   ///< pointer to plotdata

    loadStrategy loadStrategy_;

};

} // namespace voreen

#endif // VRN_PLOTENTITIESPROPERTY_H
