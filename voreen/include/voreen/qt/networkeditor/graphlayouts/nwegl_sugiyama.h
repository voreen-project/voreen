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

#ifndef VRN_NWEGL_SUGIYAMA_H
#define VRN_NWEGL_SUGIYAMA_H

#include "nwegl_base.h"

#include <vector>
#include <queue>


namespace voreen {

class VRN_QT_API NWEGL_Sugiyama : public NWEGL_Base {
public:
    NWEGL_Sugiyama();
    ~NWEGL_Sugiyama() {}

    void sort(ProcessorNetwork* network, std::vector<Processor*>* processors, QMap<Processor*,ProcessorGraphicsItem*>* itemMap);

    //set parameter
    void setSortParameter(qreal shift, bool overlap, bool median, bool portflush);
private:

    //distance between processors in same layer
    qreal shiftX_;
    //allow overlapping processors?
    bool overlap_;
    //use median positioning?
    bool median_;
    //procesors lie flush against ports?
    bool portflush_;

    //comparator for bary center sorting
    static bool pairCompare(const std::pair<double, Processor*>& firstElem, const std::pair<double, Processor*>& secondElem);

    //sorting the graphLayer based on bary center
    std::vector<std::vector<Processor*> > sortGraphLayer(std::vector<std::vector<Processor*> > graphLayer);
};

} //namespace voreen

#endif // VRN_NWEGL_FIRST_H

