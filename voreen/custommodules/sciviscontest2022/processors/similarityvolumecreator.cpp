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

#include "similarityvolumecreator.h"

#include "voreen/core/datastructures/volume/volumefactory.h"
#include "modules/ensembleanalysis/utils/utils.h"
#include "modules/plotting/datastructures/plotcell.h"
#include "modules/plotting/datastructures/plotdata.h"

#include <chrono>
#include <algorithm>
#include <random>

namespace voreen {
    SimilarityVolumeCreator::SimilarityVolumeCreator()
            : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
            , inport_(Port::INPORT, "inport", "Ensemble Input")
            , outport_(Port::OUTPORT, "outport", "RGB Volume Output")
            , eigenValueOutport_(Port::OUTPORT, "eigenvalueOutport", "Eigenvalue Output")
            , xRange_("xRange", "X-Range")
            , yRange_("yRange", "Y-Range")
            , zRange_("zRange", "Z-Range")
            , field_("field", "Field")
            , algorithm_("algorithm", "Algorithm")
            , landmarkNumber_("landmarkNumber", "Number of Landmarks", 100, 0, 10000)
            , iterations_("iterations", "Iterations")
            , member_("member", "Select member")
    {
        // add ports
        addPort(inport_);
        ON_CHANGE(inport_, SimilarityVolumeCreator, adjustToEnsemble);
        addPort(outport_);
        addPort(eigenValueOutport_);

        // properties for dimensions
        addProperty(xRange_);
        addProperty(yRange_);
        addProperty(zRange_);

        // Currently, only one member can be used
        addProperty(member_);
        // Should be changed or removed later to allow multiple fields
        addProperty(field_);

        // Properties for algorithm
        addProperty(algorithm_);
        algorithm_.addOption("Classical MDS", "Classical MDS", CLASSICAL_MDS);
        algorithm_.addOption("Landmark MDS", "Landmark MDS", LANDMARK_MDS);
        algorithm_.selectByValue(CLASSICAL_MDS);
        addProperty(landmarkNumber_);
        landmarkNumber_.setVisibleFlag(false);
        ON_CHANGE_LAMBDA(algorithm_, [this] {
            landmarkNumber_.setVisibleFlag(algorithm_.getValue() == LANDMARK_MDS);
        });

        addProperty(iterations_);
    }

    Processor* SimilarityVolumeCreator::create() const {
        return new SimilarityVolumeCreator();
    }

    SimilarityVolumeCreatorInput SimilarityVolumeCreator::prepareComputeInput() {
        PortDataPointer<EnsembleDataset> ensemble = inport_.getThreadSafeData();
        if (!ensemble) {
            throw InvalidInputException("No input", InvalidInputException::S_WARNING);
        }
        tgt::vec3 spacing = ensemble->getVolumes()[0]->getSpacing();
        float xDim = (xRange_.get().y - xRange_.get().x)/std::abs(spacing.x);
        float yDim = (yRange_.get().y - yRange_.get().x)/std::abs(spacing.y);
        float zDim = (zRange_.get().y - zRange_.get().x)/std::abs(spacing.z);
        tgt::ivec3 newDims(std::max(1, int(std::round(xDim))), std::max(1, int(std::round(yDim))), std::max(1, int(std::round(zDim))));
        VolumeFactory factory;
        // 3 channels for RGB
        std::string format = factory.getFormat("float", 3);
        std::unique_ptr<VolumeRAM> outputVolume(factory.create(format, newDims));
        outputVolume->clear(); // Sets all to zero.

        tgt::vec3 llf(xRange_.get().x, yRange_.get().x, zRange_.get().x);
        tgt::vec3 urb(xRange_.get().y, yRange_.get().y, zRange_.get().y);
        tgt::Bounds bounds(llf, urb);

        std::cout << outputVolume->getDimensions() << std::endl;

        // Clear old data.
        outport_.clear();

        return SimilarityVolumeCreatorInput {
                std::move(ensemble),
                std::move(outputVolume),
                bounds,
                field_.get(),
                algorithm_.getValue(),
                landmarkNumber_.get(),
                iterations_.get()
        };
    }

    SimilarityVolumeCreatorOutput SimilarityVolumeCreator::compute(ComputeInput input, ProgressReporter& progressReporter) const {
        auto ensemble = std::move(input.ensemble);
        const tgt::Bounds& bounds = input.bounds;
        std::string field = std::move(input.field);
        std::unique_ptr<VolumeRAM> output = std::move(input.outputVolume);
        const tgt::ivec3 newDims = output->getDimensions();

        auto start = std::chrono::high_resolution_clock::now();
        // Extract time series
        const int numVoxels = newDims.x*newDims.y*newDims.z;
        const int numTimesteps = ensemble->getMinNumTimeSteps();
        std::vector<std::vector<float>> timeseries(numVoxels, std::vector<float>(numTimesteps));
        EnsembleMember member = ensemble->getMembers()[member_.getValue()];
        for (size_t t = 0; t < ensemble->getMinNumTimeSteps(); ++t) {
            const VolumeBase* vol = member.getTimeSteps()[t].getVolume(field);
            RealWorldMapping rwm = vol->getRealWorldMapping();
            VolumeRAMRepresentationLock lock(vol);
            tgt::mat4 worldToVoxel = vol->getWorldToVoxelMatrix();
            for(int x = 0; x < newDims.x; x++) {
                for(int y = 0; y < newDims.y; y++) {
                    for(int z = 0; z < newDims.z; z++) {
                        tgt::vec3 sample = mapRange(tgt::vec3(x,y,z), tgt::vec3::zero, tgt::vec3(newDims), bounds.getLLF(), bounds.getURB());
                        tgt::svec3 sampleInVoxelSpace = worldToVoxel * sample;
                        // TODO: Only supports single channel
                        float val = rwm.normalizedToRealWorld(lock->getVoxelNormalized(sampleInVoxelSpace, 0));
                        timeseries[x*newDims.z*newDims.y+y*newDims.z+z][t] = val;
                    }
                }
            }
            progressReporter.setProgress(0.33/ensemble->getMinNumTimeSteps()*t);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
        std::cout << "Time for Time Series extraction " << time << "ms" << std::endl;

        progressReporter.setProgress(0.33f);
        // Distance matrix
        start = std::chrono::high_resolution_clock::now();
        int numPoints = timeseries.size();
        Eigen::MatrixXf PMatrix = Eigen::MatrixXf::Zero(numPoints, numPoints);
        //std::vector<float> correlations(numVoxels * numVoxels, 1.0f);
        for(int i = 0; i<numVoxels; i++) {
            for(int j = 0; j<i; j++) {
                float numerator = 0;
                float denomI = 0;
                float denomJ = 0;
                for(int t = 0; t < numTimesteps; t++) {
                    numerator += timeseries[i][t]*timeseries[j][t];
                    denomI += timeseries[i][t]*timeseries[i][t];
                    denomJ += timeseries[j][t]*timeseries[j][t];
                }
                float v = 0.5 * (1-numerator / std::sqrt(denomI * denomJ));
                PMatrix(i,j) = PMatrix(j,i) = v*v;
            }
            progressReporter.setProgress(0.33+0.33/numVoxels*i);
        }
        // Free memory
        std::vector<std::vector<float>>().swap(timeseries);
        end = std::chrono::high_resolution_clock::now();
        time = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
        std::cout << "Time for correlation calculation " << time << "ms" << std::endl;
        progressReporter.setProgress(0.66);
        // MDS
        start = std::chrono::high_resolution_clock::now();
        std::vector<float> eigenvalues;
        if(input.algorithm == CLASSICAL_MDS)
            eigenvalues = mds(PMatrix, input.iterations, *output);
        if(input.algorithm == LANDMARK_MDS)
            eigenvalues = landmarkMDS(PMatrix, input.iterations, input.k, *output);
        end = std::chrono::high_resolution_clock::now();
        time = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
        std::cout << "Time for MDS " << time << "ms" << std::endl;


        // Put coordinates in volume
        tgt::vec3 spacing = bounds.diagonal() / tgt::vec3(newDims);
        std::unique_ptr<Volume> volume(new Volume(output.release(), spacing, bounds.getLLF()));
        volume->setModality(Modality(field));

        progressReporter.setProgress(1.0f);

        return SimilarityVolumeCreatorOutput {
                std::move(volume),
                eigenvalues
        };
    }

    void SimilarityVolumeCreator::processComputeOutput(ComputeOutput output) {
        outport_.setData(output.volume.release(), true);

        //std::cout << output.eigenvalues.size() << std::endl;
        // Output eigenvalues
        PlotData* data = new PlotData(1,1);
        data->setColumnLabel(0, "Index");
        data->setColumnLabel(1, "Eigenvalue");
        for(size_t i = 0; i<output.eigenvalues.size(); i++) {
            std::vector<PlotCellValue> values;
            values.push_back(PlotCellValue(i+1));
            values.push_back(PlotCellValue(static_cast<plot_t>(output.eigenvalues[i])));
            data->insert(values);
        }
        eigenValueOutport_.setData(data, true);
    }

    void SimilarityVolumeCreator::adjustToEnsemble() {

        const EnsembleDataset* dataset = inport_.getData();
        if (!dataset)
            return;

        field_.blockCallbacks(true);
        field_.reset();
        for(const std::string& fieldName : dataset->getCommonFieldNames())
            field_.addOption(fieldName, fieldName, fieldName);
        field_.blockCallbacks(false);

        member_.reset();
        for (int i = 0; i < dataset->getMembers().size(); i++) {
            const EnsembleMember& member = dataset->getMembers()[i];
            member_.addOption(member.getName(), member.getName(), i);
        }

        tgt::vec3 llf = dataset->getBounds().getLLF();
        tgt::vec3 urb = dataset->getBounds().getURB();
        xRange_.setMinValue(llf.x);
        xRange_.setMaxValue(urb.x);
        yRange_.setMinValue(llf.y);
        yRange_.setMaxValue(urb.y);
        zRange_.setMinValue(llf.z);
        zRange_.setMaxValue(urb.z);
    }

    void SimilarityVolumeCreator::processResult(Eigen::MatrixXf result, VolumeRAM &output) const {
        float min = result.minCoeff();
        float max = result.maxCoeff();
        // Set data in volume:
        tgt::vec3 outputDimensions = output.getDimensions();
        for(int x = 0; x < outputDimensions.x; x++) {
            for(int y = 0; y < outputDimensions.y; y++) {
                for(int z = 0; z < outputDimensions.z; z++) {
                    tgt::vec3 pos(x,y,z);
                    int posResult = x*outputDimensions.y*outputDimensions.z+y*outputDimensions.z+z;
                    tgt::vec3 value(result(posResult,0),result(posResult,1),result(posResult,2));
                    value = cielabToRGB(value, min, max);
                    output.setVoxelNormalized(value.r, pos, 0);
                    output.setVoxelNormalized(value.g, pos, 1);
                    output.setVoxelNormalized(value.b, pos, 2);
                }
            }
        }
    }

    std::vector<float> SimilarityVolumeCreator::landmarkMDS(Eigen::MatrixXf &distances, int numIterations, int k,
                                                            VolumeRAM &output) const {
        using namespace Eigen;
        int numDimensions = 3;

        int numPoints = distances.cols();
        MatrixXf result(numPoints, numDimensions);

        // Select k samples
        PermutationMatrix<Dynamic, Dynamic> perm(numPoints);
        perm.setIdentity();
        std::srand(0);
        std::random_shuffle(perm.indices().data(), perm.indices().data() + perm.indices().size());
        distances = perm.transpose() * distances * perm;

        MatrixXf smallResult(k, numDimensions);
        std::vector<float> eigenvalues;
        MatrixXf eigenvectors(k, numDimensions);
        MatrixXf distanceSelection = distances.block(0,0,k,k);
        mds(distanceSelection, numIterations, smallResult, eigenvalues, eigenvectors);

        // Embed remaining points
        int N = numPoints - k;
        MatrixXf deltaN = distanceSelection.rowwise().mean().replicate(1, N);
        MatrixXf deltaX = distances.block(0, k, k, N);
        MatrixXf diff = deltaN - deltaX;
        MatrixXf furtherResult(N, numDimensions); // Result for the rest of the points
        furtherResult.setZero(); // To avoid undefined behavior
        for(int i = 0; i < numDimensions; i++) {
            if(eigenvalues[i] > 0) {
                VectorXf L = eigenvectors.col(i)/std::sqrt(eigenvalues[i]);
                furtherResult.col(i) = (0.5*L.transpose()*diff).transpose();
            }
        }
        result << smallResult, furtherResult;

        // permute back
        result = perm * result;

        processResult(result, output);

        return eigenvalues;
    }

    void SimilarityVolumeCreator::mds(Eigen::MatrixXf &distances, int numIterations, Eigen::MatrixXf &result,
                                      std::vector<float> &eigenvalues, Eigen::MatrixXf &eigenvecs) const {
        using namespace Eigen;
        int numDimensions = 3;
        float epsilon = 0.0f;

        int numPoints = distances.cols();
        MatrixXf EigSq = MatrixXf::Zero(numDimensions, numDimensions);

        // Iterative calculation of eigenvalues
        MatrixXf JMatrix = MatrixXf::Identity (numPoints, numPoints) -
                           (1.0f / numPoints) * MatrixXf::Ones (numPoints, numPoints);

        MatrixXf BMatrix = -0.5f*JMatrix*distances*JMatrix;

        //VectorXf* eigenVectors = new VectorXf[numDimensions];

        for(size_t i=0; i < numDimensions; i++) {

            VectorXf eigenVector = VectorXf::Ones(numPoints);
            float eigenValue = 0.0;

            VectorXf PrVector = VectorXf::Zero(numPoints);
            VectorXf TVector = VectorXf::Ones(numPoints);

            MatrixXf EMVector(numPoints, 1);

            for (int iter=0; iter < numIterations; iter++) {

                EMVector.col(0) = eigenVector;
                VectorXf TempVector = BMatrix * EMVector;
                eigenVector = TempVector;
                for(size_t j=0; j < i; j++)
                    eigenVector -= eigenvecs.col(j) * (eigenvecs.col(j).dot(TempVector));
                eigenValue = eigenVector.norm();
                eigenVector.normalize();
                TVector = eigenVector - PrVector;
                PrVector = eigenVector;

                if(TVector.norm() <= epsilon) {
                    break;
                }
            }

            // Don't continue calculating when eigenvalues get too small in relation to biggest eigenvalue.
            if(i >= static_cast<size_t>(numDimensions) && // We do have calculated at least as many PCs as we want to display.
               (eigenValue <= std::numeric_limits<float>::epsilon() || // Eigenvalue is de-facto zero.
                (i > 0 && eigenValue < eigenvalues[i - 1]))) // EV must not be smaller than predecessor.
                break;

            EMVector.col(0) = eigenVector;
            result.col(i) = eigenVector;
            EigSq(i, i) = std::sqrt(eigenValue);
            eigenvalues.push_back(eigenValue);
            eigenvecs.col(i) = eigenVector;
        }
        //delete [] eigenVectors;

        // Now get the resulting matrix.
        result = result * EigSq;
    }

    std::vector<float> SimilarityVolumeCreator::mds(Eigen::MatrixXf  &distances, int numIterations, VolumeRAM &output) const {
        int numDimensions = 3;
        int numPoints = distances.cols();

        Eigen::MatrixXf result(numPoints, numDimensions);
        std::vector<float> eigenvalues;
        Eigen::MatrixXf eigenvecs(numPoints, numDimensions);
        mds(distances, numIterations, result, eigenvalues, eigenvecs);

        processResult(result, output);

        return eigenvalues;

    }

    tgt::vec3 SimilarityVolumeCreator::cielabToRGB(tgt::vec3 labVector, float min, float max) const {
        labVector = labVector-min;
        float diff = max - min;
        labVector = labVector/diff*100.0f;
        labVector.y = labVector.y-50.0f;
        labVector.z = labVector.z-50.0f;
        float var_Y = ( labVector.x + 16. ) / 116.;
        float var_X = labVector.y / 500. + var_Y;
        float var_Z = var_Y - labVector.z / 200.;

        if ( pow(var_Y,3) > 0.008856 ) var_Y = pow(var_Y,3);
        else                      var_Y = ( var_Y - 16. / 116. ) / 7.787;
        if ( pow(var_X,3) > 0.008856 ) var_X = pow(var_X,3);
        else                      var_X = ( var_X - 16. / 116. ) / 7.787;
        if ( pow(var_Z,3) > 0.008856 ) var_Z = pow(var_Z,3);
        else                      var_Z = ( var_Z - 16. / 116. ) / 7.787;

        float X = 95.047 * var_X ;    //ref_X =  95.047     Observer= 2°, Illuminant= D65
        float Y = 100.000 * var_Y  ;   //ref_Y = 100.000
        float Z = 108.883 * var_Z ;    //ref_Z = 108.883


        var_X = X / 100. ;       //X from 0 to  95.047      (Observer = 2°, Illuminant = D65)
        var_Y = Y / 100. ;       //Y from 0 to 100.000
        var_Z = Z / 100. ;      //Z from 0 to 108.883

        //std::cout << var_X << " " << var_Y << " " << var_Z << std::endl;

        float var_R = var_X *  3.2406 + var_Y * -1.5372 + var_Z * -0.4986;
        float var_G = var_X * -0.9689 + var_Y *  1.8758 + var_Z *  0.0415;
        float var_B = var_X *  0.0557 + var_Y * -0.2040 + var_Z *  1.0570;

        if ( var_R > 0.0031308 ) var_R = 1.055 * pow(var_R , ( 1 / 2.4 ))  - 0.055;
        else                     var_R = 12.92 * var_R;
        if ( var_G > 0.0031308 ) var_G = 1.055 * pow(var_G , ( 1 / 2.4 ) )  - 0.055;
        else                     var_G = 12.92 * var_G;
        if ( var_B > 0.0031308 ) var_B = 1.055 * pow( var_B , ( 1 / 2.4 ) ) - 0.055;
        else                     var_B = 12.92 * var_B;

        tgt::vec3 rgb;
        rgb.r = std::min(std::max(var_R * 255., 0.0), 255.0);
        rgb.g = std::min(std::max(var_G * 255., 0.0), 255.0);
        rgb.b = std::min(std::max(var_B * 255., 0.0), 255.0);
        return rgb;
    }

}
