/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_ENSEMBLEFILTER_H
#define VRN_ENSEMBLEFILTER_H

#include "voreen/core/processors/processor.h"

#include "../ports/ensembledatasetport.h"

namespace voreen {

class Filter {
public:
    virtual ~Filter() {}
    virtual Property& getProperty() = 0;
    virtual EnsembleDataset* applyFilter(const EnsembleDataset* ensemble) = 0;
    virtual void adjustToEnsemble(const EnsembleDataset* ensemble) = 0;
};

/**
 * Base class for all processors filtering an ensemble dataset.
 */
class VRN_CORE_API EnsembleFilter : public Processor {
public:
    EnsembleFilter();
    virtual ~EnsembleFilter();

    virtual Processor* create() const;
    virtual std::string getClassName() const { return "EnsembleFilter";        }
    virtual std::string getCategory() const  { return "Filter";                }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL; }
    virtual bool isUtility() const           { return true;                    }

protected:

    void addFilter(Filter* filter);

    void process();
    void invalidate(int inv = 1);

    void applyFilter();
    void adjustToEnsemble();

    /// Inport for the ensemble data structure.
    EnsembleDatasetPort ensembleInport_;

    /// The ensemble data
    EnsembleDatasetPort ensembleOutport_;

    /// Filter list
    std::vector<Filter*> filters_;

    /// Determines if process needs to be executed.
    bool needsProcess_;
};

} // namespace

#endif