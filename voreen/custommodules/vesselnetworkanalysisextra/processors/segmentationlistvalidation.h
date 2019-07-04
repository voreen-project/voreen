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

#ifndef VRN_SEGMENTATIOLISTNVALIDATION_H
#define VRN_SEGMENTATIOLISTNVALIDATION_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/filedialogproperty.h"

namespace voreen {

/**
 * Measures the quality of a list of volume segmentations by comparing it against a reference segmentation.
 *
 * Several similarity measures are computed, such as Jaccard and Dice indices. The results are
 * displayed by read-only properties.
 *
 * Note: This is heavily based on the original SegmentationValidation, which only compares one segmentation at a time
 *
 * !!!!!
 * Note: This processor currently includes a hack that assumes that the threshold used for segmentation is stored as
 *       offset in the RWM!
 * !!!!!
 */
class VRN_CORE_API SegmentationListValidation : public Processor {
public:
    SegmentationListValidation();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "SegmentationListValidation"; }
    virtual std::string getCategory() const  { return "Utility"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL; }
    virtual bool isUtility() const           { return true; }
    virtual bool isEndProcessor() const      { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Measures the quality of a volume segmentation list by comparing it against a reference segmentation. Several similarity measures are computed.");
    }

    virtual void process();

private:
    void forceUpdate();
    void resetOutput();

    void computeMetrics();
    size_t getNumCommonVoxels(const VolumeRAM* volumeA, const VolumeRAM* volumeB);

    VolumeListPort inportSegmentations_;
    VolumePort inportReference_;

    FileDialogProperty exportFilePath_;

    ButtonProperty computeButton_;

    bool forceUpdate_;

private:
    static const std::string loggerCat_; ///< category used in logging
};


} // namespace

#endif // VRN_SEGMENTATIONVALIDATION_H
