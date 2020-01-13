/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "linkevaluatorplotselection.h"
#include "../plotselectionproperty.h"
#include "../../datastructures/plotselection.h"

namespace voreen {

void LinkEvaluatorPlotSelection::eval(Property* src, Property* dst) {
    std::vector<PlotSelectionEntry> targetNewVec = static_cast<PlotSelectionProperty*>(src)->get();
    std::vector<PlotSelectionEntry> targetOldVec = static_cast<PlotSelectionProperty*>(dst)->get();
    std::vector<PlotSelectionEntry>::iterator targetNewIt;
    std::vector<PlotSelectionEntry>::iterator targetOldIt;
    for (targetNewIt = targetNewVec.begin(); targetNewIt < targetNewVec.end(); ++targetNewIt) {
        // check if *targetNewIt is a new entry
        bool newEntry = true;
        for (targetOldIt = targetOldVec.begin(); targetOldIt < targetOldVec.end(); ++targetOldIt) {
            if (targetOldIt->selection_ == targetNewIt->selection_) {
                newEntry = false;
                break;
            }
        }
        if (newEntry) {
            targetNewIt->highlight_ = false;
            targetNewIt->zoomTo_ = false;
            targetNewIt->renderLabel_ = false;
        }
        else {
            targetNewIt->highlight_ = targetOldIt->highlight_;
            targetNewIt->zoomTo_ = targetOldIt->zoomTo_;
            targetNewIt->renderLabel_ = targetOldIt->renderLabel_;
        }
    }
    static_cast<PlotSelectionProperty*>(dst)->set(targetNewVec);
}

bool LinkEvaluatorPlotSelection::arePropertiesLinkable(const Property* p1, const Property* p2) const {
    tgtAssert(p1, "null pointer");
    tgtAssert(p2, "null pointer");

    return (dynamic_cast<const PlotSelectionProperty*>(p1) && dynamic_cast<const PlotSelectionProperty*>(p2));
}

} // namespace
