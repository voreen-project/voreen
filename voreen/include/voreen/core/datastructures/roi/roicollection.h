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

#ifndef VRN_ROICOLLECTION_H
#define VRN_ROICOLLECTION_H

#include "voreen/core/datastructures/roi/roibase.h"

namespace voreen {

class ROICollection: public Serializable {
    public:
    ~ROICollection();

    size_t getNumROIs() const { return rois_.size(); }
    ROIBase* getROI(size_t i) { return rois_[i]; }
    const ROIBase* getROI(size_t i) const { return rois_[i]; }
    void addROI(ROIBase* roi) { rois_.push_back(roi); }
    void removeROI(ROIBase* roi);
    void deregisterROI(ROIBase* roi);
    void clear();
    virtual bool moveControlPoint(const ControlPoint* cp, tgt::vec3 to);

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    void save(const std::string& filename) const;
    void load(const std::string& filename);

    private:
    bool removeROIrecursive(ROIBase* p, ROIBase* roi);
    std::vector<ROIBase*> rois_;

    static const std::string loggerCat_;
};

} //namespace

#endif
