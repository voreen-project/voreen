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

#ifndef VRN_MUTUALINFORMATIONREGISTRATION_H
#define VRN_MUTUALINFORMATIONREGISTRATION_H

#include "itkprocessor.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/properties/optionproperty.h"

#include <string>
#include <boost/thread.hpp>

#include "itkImage.h"
#include "itkVersorRigid3DTransform.h"
#include "itkVersorRigid3DTransformOptimizer.h"

namespace voreen {

class Volume;
class MutualInformationRegistration;

//  The following section of code implements an observer
//  that will monitor the evolution of the registration process.
class MutualInformationRegistrationObserver : public itk::Command
{
    public:
        typedef  MutualInformationRegistrationObserver   Self;
        typedef  itk::Command             Superclass;
        typedef  itk::SmartPointer<Self>  Pointer;
        itkNewMacro( Self );
    protected:
        MutualInformationRegistration* processor_;
        MutualInformationRegistrationObserver() : processor_(0) {};
    public:
        typedef   itk::VersorRigid3DTransformOptimizer OptimizerType;

        void setProcessor(MutualInformationRegistration* processor) { processor_ = processor; }

        void Execute(itk::Object *caller, const itk::EventObject & event);
        //{
            //Execute( (const itk::Object *)caller, event);
        //}

        void Execute(const itk::Object * object, const itk::EventObject & event) {
        }
};

class MutualInformationRegistration : public VolumeProcessor {
    typedef itk::VersorRigid3DTransform<double> TransformType;
    typedef TransformType::ParametersType ParametersType;
public:
    friend class MutualInformationRegistrationObserver;

    MutualInformationRegistration();
    virtual Processor* create() const;

    virtual std::string getCategory() const   { return "Volume Processing"; }
    virtual std::string getClassName() const  { return "MutualInformationRegistration";  }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }

protected:
    virtual void setDescriptions() {
        setDescription("Computes a coregistration matrix for two volumes using the multi-resolution framework implemented in ITK.");
    }

    virtual void deinitialize();

    virtual void process();

    void stopRegistration();
    void updateRegistration();
    template<class fixedType>
    void typeSwitch(const VolumeBase* fixedVolume, const VolumeBase* movingVolume);
    template<class fixedType, class movingType>
    void performRegistration(const VolumeBase* fixedVolume, const VolumeBase* movingVolume);

    void initializeRegistration();
    template<class fixedType>
    void typeSwitchInit(const VolumeBase* fixedVolume, const VolumeBase* movingVolume);
    template<class fixedType, class movingType>
    void performInitialization(const VolumeBase* fixedVolume, const VolumeBase* movingVolume);

    void onMetricChange();

    tgt::mat4 calculateVoreenTrafo(const tgt::mat4& itkMatrix);
    tgt::mat4 calculateITKTrafo();
    tgt::mat4 itkToVoreen(TransformType::Pointer t);

    const VolumeBase* getFixedVolume() const;
    const VolumeBase* getMovingVolume() const;

private:
    VolumePort fixedVolumeInport_;
    VolumePort movingVolumeInport_;

    FloatMat4Property transformationMatrix_;
    IntProperty numLevels_;
    IntProperty numIterations_;
    FloatProperty maxStepLength_;
    FloatProperty minStepLength_;
    FloatProperty rotScale_;
    FloatProperty relaxationFactor_;

    // Metric properties:
    StringOptionProperty metric_;
    IntProperty numHistogramBins_;
    IntProperty numSamples_;
    BoolProperty explicitPDF_;
    FloatProperty lambda_;

    ButtonProperty updateButton_;
    ButtonProperty stopButton_;
    ButtonProperty initializeButton_;

    boost::thread workerThread_;
    tgt::mat4 tempResult_;
    bool tempResultUpdated_;

    static const std::string loggerCat_; ///< category used in logging
};

}   //namespace

#endif
