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

#include "itkprocessor.h"

#include "itkProcessObject.h"

namespace voreen {

const std::string ITKProcessor::loggerCat_("voreen.ITKProcessor");

ITKProcessor::ITKProcessor()
    : CachingVolumeProcessor()
{
    redrawCommand_ = RedrawCommandType::New();
    redrawCommand_->SetCallbackFunction(this, &ITKProcessor::ProcessEvent);
    redrawCommand_->SetCallbackFunction(this, &ITKProcessor::ConstProcessEvent);
}

void ITKProcessor::ProcessEvent(itk::Object* caller, const itk::EventObject& event) {
    if(typeid(itk::ProgressEvent) == typeid(event))
    {
        ::itk::ProcessObject::Pointer process = dynamic_cast<itk::ProcessObject*>(caller);
        setProgress(process->GetProgress());
    }
}

void ITKProcessor::ConstProcessEvent(const itk::Object* caller, const itk::EventObject& event) {
    if(typeid(itk::ProgressEvent) == typeid(event))
    {
        itk::ProcessObject::ConstPointer process = dynamic_cast<const itk::ProcessObject*>(caller);
        setProgress(process->GetProgress());
    }
}

void ITKProcessor::observe(itk::Object* caller) {
      caller->AddObserver(itk::ProgressEvent(), redrawCommand_.GetPointer());
}

}   // namespace
