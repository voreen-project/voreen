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

#ifndef VRN_SEGMENTATIONVALIDATION_H
#define VRN_SEGMENTATIONVALIDATION_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/boolproperty.h"

namespace voreen {

/**
 * Measures the quality of a volume segmentation by comparing it against a reference segmentation.
 *
 * Several similarity measures are computed, such as Jaccard and Dice indices. The results are
 * displayed by read-only properties.
 */
class VRN_CORE_API SegmentationValidation : public Processor {
public:
    SegmentationValidation();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "SegmentationValidation"; }
    virtual std::string getCategory() const  { return "Utility"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE; }
    virtual bool isUtility() const           { return true; }
    virtual bool isEndProcessor() const      { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Measures the quality of a volume segmentation by comparing it against a reference segmentation. Several similarity measures are computed. The results are displayed by read-only properties.");
    }
    
    virtual void beforeProcess();
    
    virtual void process();

private:
    void forceUpdate();
    void resetOutput();

    void computeMetrics();
    size_t getNumCommonVoxels(const VolumeRAM* volumeA, const VolumeRAM* volumeB);

    /// exports quantification data to a CSV file
    void exportData();

    VolumePort inportSegmentation_;
    VolumePort inportReference_;

    ButtonProperty computeButton_;
    IntProperty sizeSegmentation_;
    IntProperty sizeReference_;
    IntProperty truePositive_;
    IntProperty falsePositive_;
    IntProperty falseNegative_;
    IntProperty trueNegative_;
    FloatProperty accuracy_;
    FloatProperty accuracyForeground_;
    FloatProperty accuracyBackground_;
    FloatProperty jaccardIndex_;
    FloatProperty diceIndex_;
    FloatProperty sensitivity_;
    FloatProperty specificity_;

    // export data properties
    FileDialogProperty exportFile_;
    ButtonProperty exportButton_;

    // auto-computation for automatic evaluation
    BoolProperty autoCompute_;
    BoolProperty autoExport_;

    bool forceUpdate_;

private:
    static const std::string loggerCat_; ///< category used in logging
};


} // namespace

#endif // VRN_SEGMENTATIONVALIDATION_H
