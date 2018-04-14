/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#ifndef VRN_ITKPROCESSOR_H
#define VRN_ITKPROCESSOR_H

#include "voreen/core/processors/volumeprocessor.h"
#include "itkCommand.h"

namespace voreen {

class Volume;

class ITKProcessor : public CachingVolumeProcessor {
public:
    ITKProcessor();

    virtual bool usesExpensiveComputation() const { return true; }

    //Stuff to get progress from ITK:

    /// Manage a Progress event
    void ProcessEvent(itk::Object* caller, const itk::EventObject& event);
    void ConstProcessEvent(const itk::Object* caller, const itk::EventObject& event);

    /// Observe an ITK object (usually a filter)
    void observe(itk::Object* caller);

private:
    /// Command Class invoked for progress updates
    typedef itk::MemberCommand<ITKProcessor> RedrawCommandType;

    RedrawCommandType::Pointer redrawCommand_;

    static const std::string loggerCat_; ///< category used in logging
};

}   //namespace

#endif
