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

#include "lambdacicriterion.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorgradient.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"


#include "Eigen/Eigenvalues"

namespace voreen {

const std::string LambdaCiCriterion::loggerCat_("voreen.flowsimulation.lambdacicriterion");

LambdaCiCriterion::LambdaCiCriterion()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inputVolume_(Port::INPORT, "lambdacicriterion.inputVolume", "Volume Input")
    , outputVolume_(Port::OUTPORT, "lambdacicriterion.outputVolume", "Volume Output")
{
    addPort(inputVolume_);
    inputVolume_.addCondition(new PortConditionVolumeChannelCount(3));
    addPort(outputVolume_);
}

Processor* LambdaCiCriterion::create() const {
    return new LambdaCiCriterion();
}

LambdaCiCriterionInput LambdaCiCriterion::prepareComputeInput() {

    auto inputVolume = inputVolume_.getThreadSafeData();
    if(!inputVolume) {
        throw InvalidInputException("No input", InvalidInputException::S_ERROR);
    }

    std::unique_ptr<VolumeRAM_Float> outputVolume(new VolumeRAM_Float(inputVolume->getDimensions()));

    return LambdaCiCriterionInput{
        std::move(inputVolume),
        std::move(outputVolume)
    };
}
/*
void getEigenValues(float& l1, float& l2, float& l3) {
    // Set up characteristic equation:   det( A - lambda I ) = 0
    //    as a cubic in lambda:  a.lambda^3 + b.lambda^2 + c.lambda + d = 0
    std::complex<double> a = -1.0;                 // -1
    std::complex<double> b = xx + yy + zz;         // trace
    std::complex<double> c = yz * yz - yy * zz     // -sum of diagonal minors
                             + xz * xz - zz * xx
                             + xy * xy - xx * yy;
    std::complex<double> d = xx*yy*zz + 2*xy*yz*xz - xx*yz*yz - yy*xz*xz - zz*xy*xy;

    // Solve cubic by Cardano's method (easier in complex numbers!)
    std::complex<double> p = ( b * b - 3.0 * a * c ) / ( 9.0 * a * a );
    std::complex<double> q = ( 9.0 * a * b * c - 27.0 * a * a * d - 2.0 * b * b * b ) / ( 54.0 * a * a * a );
    std::complex<double> delta = q * q - p * p * p;
    std::complex<double> g1 = pow( q + sqrt( delta ), 1.0 / 3.0 );     // warning: exponents and sqrt of complex
    std::complex<double> g2 = pow( q - sqrt( delta ), 1.0 / 3.0 );
    std::complex<double> offset = -b / ( 3.0 * a );
    std::complex<double> omega = std::complex<double>( -0.5, 0.5 * sqrt( 3.0 ) );     // complex cube root of unity
    std::complex<double> omega2 = omega * omega;

    std::complex<double> cl1 = g1          + g2          + offset;
    std::complex<double> cl2 = g1 * omega  + g2 * omega2 + offset;
    std::complex<double> cl3 = g1 * omega2 + g2 * omega  + offset;

    //tgtAssert(std::abs(cl1.imag()) <= std::abs(cl1.real()), "Invalid eigenvalue for symmat");
    //tgtAssert(std::abs(cl2.imag()) <= std::abs(cl2.real()), "Invalid eigenvalue for symmat");
    //tgtAssert(std::abs(cl3.imag()) <= std::abs(cl3.real()), "Invalid eigenvalue for symmat");

    l1=cl1.real();
    l2=cl2.real();
    l3=cl3.real();
}
*/
LambdaCiCriterionOutput LambdaCiCriterion::compute(LambdaCiCriterionInput input, ProgressReporter& progressReporter) const {

    auto inputVolume = std::move(input.inputVolume);
    VolumeRAMRepresentationLock inputVolumeData(inputVolume);
    auto outputVolume = std::move(input.outputVolume);
    outputVolume->clear();

    tgt::svec3 dimensions = inputVolume->getDimensions();
    tgt::vec3 spacing = inputVolume->getSpacing();

    auto jacobi = Eigen::Matrix3f();

    // Very ugly, but allows to easily switch to different gradient function.
    //VolumeOperatorGradient operatorGradient;
    auto gradientFunc = std::bind(
            static_cast<tgt::vec3(&)(const VolumeRAM*, const tgt::vec3&, const tgt::svec3&, size_t)>(
                    VolumeOperatorGradient::calcGradientCentralDifferences
                    ),
            std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4
    );

    ThreadedTaskProgressReporter progress(progressReporter, dimensions.z);

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
    for (size_t z = 0; z < dimensions.z; z++) {
        for (size_t y = 0; y < dimensions.y; y++) {
            for (size_t x = 0; x < dimensions.x; x++) {
                tgt::svec3 pos(x, y, z);
                for(size_t i = 0; i < 3; i++) {
                    tgt::vec3 gradient = gradientFunc(*inputVolumeData, spacing, pos, i);
                    for(size_t j = 0; j < 3; j++) {
                        jacobi(i, j) = gradient[j];
                    }
                }

                const float epsilon = 0.00001f;

                auto eigenvalues = jacobi.eigenvalues();

                std::vector<int> real;
                std::vector<int> imag;
                for(size_t i=0; i<3; i++) {
                    if(std::abs(eigenvalues(i).imag()) > epsilon) {
                        imag.push_back(i);
                    }
                    else {
                        real.push_back(i);
                    }
                }

                if(real.size() != 1) {
                    //std::cout << "Only one real allowed!" << std::endl;
                    continue;
                }

                if(std::abs(eigenvalues(imag[0]).real() - eigenvalues(imag[1]).real()) > epsilon) {
                    //std::cout << "Real parts of imaginary eigenvaules differ!" << std::endl;
                    continue;
                }

                if(std::abs(eigenvalues(imag[0]).imag() + eigenvalues(imag[1]).imag()) > epsilon) {
                    //std::cout << "Image parts of imaginary eigenvaules differ!" << std::endl;
                    continue;
                }

                outputVolume->voxel(pos) = 1.0f / std::abs(eigenvalues(imag[0]).imag());
                //std::cout << "Found!" << std::endl;

                /*
                bool invalid = false;
                
                boost::optional<int> cr = boost::none;
                boost::optional<int> ci = boost::none;

                for(size_t i=0; i<3 && !invalid; i++) {

                    float imag = eigenvalues[i].imag();
                    if(std::abs(imag) < epsilon) { // Found real valued eigenvalue.
                        if (cr) {
                            invalid = true; // Already found one, so criterion is not met.
                            break;
                        }
                        else {
                            cr = i; // Set real eigenvalue idx.
                        }
                    }
                    else if(std::abs(imag) > epsilon) { // Found complex eigenvalue.
                        if(ci) {
                            // Second complex eigenvalue, check if it's complex conjugated of first one.
                            float real = eigenvalues[i].real(); // Real values must be (mostly) equal, so pick any as reference.
                            if(std::abs(eigenvalues[*ci].real() - real) > epsilon ||
                               std::abs(eigenvalues[*ci].imag() + imag) > epsilon)
                            {
                                invalid = true;
                                break;
                            }
                        }
                        else {
                            ci = i; // First complex eigenvalue, just set it.
                        }
                    }
                }

                if(!invalid) {
                    outputVolume->voxel(pos) = 1.0f / std::abs(eigenvalues[*ci].imag());
                }

                */
            }
        }

        progress.reportStepDone();
    }

    std::unique_ptr<Volume> volume(new Volume(outputVolume.release(), inputVolume));
    volume->setMetaDataValue<StringMetaData>("name", "lambdaCi");

    progressReporter.setProgress(1.0f);

    return LambdaCiCriterionOutput{
        std::move(volume)
    };
}

void LambdaCiCriterion::processComputeOutput(LambdaCiCriterionOutput output) {
    outputVolume_.setData(output.volume.release());
}

} // namespace voreen
