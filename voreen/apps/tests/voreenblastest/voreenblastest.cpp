/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "voreen/core/utils/voreenblas/voreenblascpu.h"
#include "modules/opencl/utils/voreenblascl.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/commandlineparser.h"
#include "modules/opencl/openclmodule.h"

#include "tgt/logmanager.h"
#include "tgt/init.h"

using namespace voreen;

inline int volumeCoordsToIndex(int x, int y, int z, const tgt::ivec3& dim) {
    return z*dim.y*dim.x + y*dim.x + x;
}

bool testSAXPY();
bool testSDOT();
bool testSNRM2();
bool testSSpMVEll();
bool testHSpMVEll();
bool testSSpInnerProductEll();
bool testSSpConjGradEll();
bool testHSpConjGradEll();

void randomizeEll(EllpackMatrix<float>& mat);
void randomizeEll(EllpackMatrix<int16_t>& mat);
void randomizeEllPositiveDefinite(EllpackMatrix<float>& mat);
void randomizeEllPositiveDefinite(EllpackMatrix<int16_t>& mat);
bool cmpVectors(size_t vecSize, float* vecx, float* vecy, float relThresh = 1e6);

const std::string loggerCat_ = "VoreenBlasTest";

VoreenBlasCPU voreenBlasCPU;
VoreenBlasCL voreenBlasCL;

const unsigned int DEFAULT_MATRIX_SEED = static_cast<unsigned int>(time(0));
int MAX_VECTOR_SIZE = 20;
int MAX_MATRIX_SIZE = 16;
float ERROR_THRESH;
unsigned int MATRIX_SEED;

int main(int argc, char* argv[]) {

    VoreenApplication app("voreenblastest", "VoreenBlasTest", "", argc, argv, VoreenApplication::APP_ALL);
    bool quickMode = false;
    app.getCommandLineParser()->addFlagOption("quickmode", quickMode, CommandLineParser::MainOption, "Execute tests with smaller data sizes");
    app.getCommandLineParser()->addOption("errorThreshold", ERROR_THRESH, CommandLineParser::MainOption, "Deviation threshold whose exceeding defines a test failure", 1e-6f, "1e-6f");
    app.getCommandLineParser()->addOption("matrixSeed", MATRIX_SEED, CommandLineParser::MainOption, "Seed to be used for randomized matrix generation", DEFAULT_MATRIX_SEED, "randomized");
    app.initialize();

    // Set seed.
    LINFO("Using seed " << MATRIX_SEED << " for randomized tests");
    srand(MATRIX_SEED);

    // Test data size.
    if (quickMode) {
        MAX_VECTOR_SIZE = 20;
        MAX_MATRIX_SIZE = 15;
        LINFO("Running in quick mode: MAX_VECTOR_SIZE=2^" << MAX_VECTOR_SIZE << ", MAX_MATRIX_SIZE=2^" << MAX_MATRIX_SIZE);
    }
    else {
        MAX_VECTOR_SIZE = 24;
        MAX_MATRIX_SIZE = 20;
        LINFO("Running in normal mode: MAX_VECTOR_SIZE=2^" << MAX_VECTOR_SIZE << ", MAX_MATRIX_SIZE=2^" << MAX_MATRIX_SIZE);
    }

    // Disable GL sharing in OpenCL module, since we do not have a OpenGL context.
    if (OpenCLModule::getInstance()) {
        OpenCLModule::getInstance()->setGLSharing(false);
    }
    else {
        LERROR("OpenCLModule not instantiated.");
        return 1;
    }

    // Initialize OpenCL.
    try {
        OpenCLModule::getInstance()->initializeGL();
    }
    catch (VoreenException& e) {
        LERROR("Failed to initialize OpenCL: " << e.what());
        return 1;
    }

    // Initialize VoreenBlasCL.
    try {
        voreenBlasCL.initialize();
    }
    catch (VoreenException& e) {
        LERROR("Failed to initialize VoreenBlasCL: " << e.what());
        return 1;
    }

    // run tests
    int numTests = 0;
    int success = 0;
    clock_t start = clock();

    success += testSDOT() ? 1 : 0;
    numTests++;

    success += testSNRM2() ? 1 : 0;
    numTests++;

    success += testSAXPY() ? 1 : 0;
    numTests++;

    success += testSSpMVEll() ? 1 : 0;
    numTests++;

    success += testHSpMVEll() ? 1 : 0;
    numTests++;

    success += testSSpInnerProductEll() ? 1 : 0;
    numTests++;

    success += testSSpConjGradEll() ? 1 : 0;
    numTests++;

    success += testHSpConjGradEll() ? 1 : 0;
    numTests++;

    clock_t end = clock();
    float duration = static_cast<float>(end-start)/CLOCKS_PER_SEC;

    // Deinitialize VoreenBlasCL.
    try {
        voreenBlasCL.deinitialize();
    }
    catch (VoreenException& e) {
        LERROR("Failed to deinitialize VoreenBlasCL: " << e.what());
        return 1;
    }

    // Deinitialize OpenCL.
    try {
        OpenCLModule::getInstance()->deinitializeGL();
    }
    catch (VoreenException& e) {
        LERROR("Failed to deinitialize OpenCL: " << e.what());
        return 1;
    }

    if (success == numTests) {
        LINFO("SUCCESS: " << success << "/" << numTests << " sub tests passed.");
        LINFO("Test Duration: " << duration << " sec");
        app.deinitialize();
        return 0;
    }
    else {
        LWARNING("FAILED: " << success << "/" << numTests << " sub tests passed. See test results above for details.");
        LINFO("Test Duration: " << duration << " sec");
        app.deinitialize();
        return 1;
    }
}

bool testSAXPY() {

    int maxDim = (1<<MAX_VECTOR_SIZE) + 17;
    LINFO("sAXPY Randomized Tests:");

    // initialize test vectors
    float* x = new float[maxDim];
    float* y = new float[maxDim];
    float* resultCL = new float[maxDim];
    float* resultCPU = new float[maxDim];
    float* diffBuff = new float[maxDim];
    for (int i=0; i<maxDim; i++) {
        x[i] = ((float)rand() / RAND_MAX) - 0.5f;
        y[i] = ((float)rand() / RAND_MAX) * 2.f - 0.6f;
    }

    int numTests = 0;
    int numSucc = 0;

    // run tests for varying vector sizes
    for (int i=2; i<=MAX_VECTOR_SIZE; ++i) {

        size_t n = (1<<i);
        float alpha = 1.f;

        // compare results computed on GPU and CPU
        voreenBlasCL.sAXPY(n, x, y, alpha, resultCL);
        voreenBlasCPU.sAXPY(n, x, y, alpha, resultCPU);
        voreenBlasCPU.sAXPY(n, resultCL, resultCPU, -1.f, diffBuff);
        float norm = voreenBlasCPU.sNRM2(n, diffBuff);
        if (norm < ERROR_THRESH*n) {
            LINFO("n=" << n << " alpha=" << alpha <<  " : " << "passed (" << norm << ")");
            numSucc++;
        }
        else
            LWARNING("n=" << n << " alpha=" << alpha <<  " : " << "FAILED  (" << norm << ")");
        numTests++;

        // repeat test for non-power-of-two vector sizes and varying alphas
        n = (1<<i) - 1;
        alpha = 0.5f;
        voreenBlasCL.sAXPY(n, x, y, alpha, resultCL);
        voreenBlasCPU.sAXPY(n, x, y, alpha, resultCPU);
        voreenBlasCPU.sAXPY(n, resultCL, resultCPU, -1.f, diffBuff);
        norm = voreenBlasCPU.sNRM2(n, diffBuff);
        if (norm < ERROR_THRESH*n) {
            LINFO("n=" << n << " alpha=" << alpha <<  " : " << "passed (" << norm << ")");
            numSucc++;
        }
        else
            LWARNING("n=" << n << " alpha=" << alpha <<  " : " << "FAILED  (" << norm << ")");
        numTests++;

        n = (1<<i) + 17;
        alpha = -2.5f;
        voreenBlasCL.sAXPY(n, x, y, alpha, resultCL);
        voreenBlasCPU.sAXPY(n, x, y, alpha, resultCPU);
        voreenBlasCPU.sAXPY(n, resultCL, resultCPU, -1.f, diffBuff);
        norm = voreenBlasCPU.sNRM2(n, diffBuff);
        if (norm < ERROR_THRESH*n) {
            LINFO("n=" << n << " alpha=" << alpha <<  " : " << "passed (" << norm << ")");
            numSucc++;
        }
        else
            LWARNING("n=" << n << " alpha=" << alpha <<  " : " << "FAILED  (" << norm << ")");
        numTests++;

    }

    delete[] x;
    delete[] y;
    delete[] resultCL;
    delete[] resultCPU;
    delete[] diffBuff;

    if (numSucc == numTests) {
        LINFO("sAXPY Randomized Tests (passed/overall): " << numSucc << "/" << numTests);
        return true;
    }
    else {
        LWARNING("sAXPY Randomized Tests (passed/overall): " << numSucc << "/" << numTests);
        return false;
    }
}

bool testSDOT() {

    int maxDim = (1<<MAX_VECTOR_SIZE) + 17;

    LINFO("sDOT Randomized Tests:");

    // initialize test vectors
    float* x = new float[maxDim];
    float* y = new float[maxDim];
    for (int i=0; i<maxDim; i++) {
        x[i] = ((float)rand() / RAND_MAX) - 0.5f;
        y[i] = ((float)rand() / RAND_MAX) * 2.f - 0.6f;
    }

    int numTests = 0;
    int numSucc = 0;

    // run tests for varying vector sizes
    for (int i=2; i<=MAX_VECTOR_SIZE; ++i) {

        size_t n = (1<<i);

        // compare results computed on GPU and CPU
        float resultCL = voreenBlasCL.sDOT(n, x, y);
        float resultCPU = voreenBlasCPU.sDOT(n, x, y);
        float relError = std::abs((resultCL / resultCPU) - 1.f);
        if (relError < std::min(ERROR_THRESH*n, 0.01f)) {
            LINFO("n=" << n << ":  Result GPU: " << resultCL << ", Result CPU: " << resultCPU << ", Relative Error: " << relError);
            numSucc++;
        }
        else
            LWARNING("n=" << n << ":  Result GPU: " << resultCL << ", Result CPU: " << resultCPU << ", Relative Error: " << relError << " (FAILED)");
        numTests++;

        // repeat test for non-power-of-two vector sizes
        n = (1<<i) - 1;
        resultCL = voreenBlasCL.sDOT(n, x, y);
        resultCPU = voreenBlasCPU.sDOT(n, x, y);
        relError = std::abs((resultCL / resultCPU) - 1.f);
        if (relError < ERROR_THRESH*n) {
            LINFO("n=" << n << ":  Result GPU: " << resultCL << ", Result CPU: " << resultCPU << ", Relative Error: " << relError);
            numSucc++;
        }
        else
            LWARNING("n=" << n << ":  Result GPU: " << resultCL << ", Result CPU: " << resultCPU << ", Relative Error: " << relError << " (FAILED)");
        numTests++;

        n = (1<<i) + 17;
        resultCL = voreenBlasCL.sDOT(n, x, y);
        resultCPU = voreenBlasCPU.sDOT(n, x, y);
        relError = std::abs((resultCL / resultCPU) - 1.f);
        if (relError < std::min(ERROR_THRESH*n, 0.01f)) {
            LINFO("n=" << n << ":  Result GPU: " << resultCL << ", Result CPU: " << resultCPU << ", Relative Error: " << relError);
            numSucc++;
        }
        else
            LWARNING("n=" << n << ":  Result GPU: " << resultCL << ", Result CPU: " << resultCPU << ", Relative Error: " << relError << " (FAILED)");
        numTests++;

    }

    delete[] x;
    delete[] y;

    if (numSucc == numTests) {
        LINFO("sDOT Randomized Tests (passed/overall): " << numSucc << "/" << numTests);
        return true;
    }
    else {
        LWARNING("sDOT Randomized Tests (passed/overall): " << numSucc << "/" << numTests);
        return false;
    }
}

bool testSNRM2() {

    int maxDim = (1<<MAX_VECTOR_SIZE) + 17;
    LINFO("sNRM2 Randomized Tests:");

    // initialize test vectors
    float* x = new float[maxDim];
    for (int i=0; i<maxDim; i++) {
        x[i] = ((float)rand() / RAND_MAX) - 0.5f;
    }

    int numTests = 0;
    int numSucc = 0;

    // run tests for varying vector sizes
    for (int i=2; i<=MAX_VECTOR_SIZE; ++i) {

        size_t n = (1<<i);

        // compare results computed on GPU and CPU
        float resultCL = voreenBlasCL.sNRM2(n, x);
        float resultCPU = voreenBlasCPU.sNRM2(n, x);
        float relError = std::abs((resultCL / resultCPU) - 1.f);
        if (relError < std::min(ERROR_THRESH*n, 0.1f)) {
            LINFO("n=" << n << ":  Result GPU: " << resultCL << ", Result CPU: " << resultCPU << ", Relative Error: " << relError);
            numSucc++;
        }
        else
            LWARNING("n=" << n << ":  Result GPU: " << resultCL << ", Result CPU: " << resultCPU << ", Relative Error: " << relError << " (FAILED)");
        numTests++;

        // repeat test for non-power-of-two vector sizes
        n = (1<<i) - 1;
        resultCL = voreenBlasCL.sNRM2(n, x);
        resultCPU = voreenBlasCPU.sNRM2(n, x);
        relError = std::abs((resultCL / resultCPU) - 1.f);
        if (relError < std::min(ERROR_THRESH*n, 0.1f)) {
            LINFO("n=" << n << ":  Result GPU: " << resultCL << ", Result CPU: " << resultCPU << ", Relative Error: " << relError);
            numSucc++;
        }
        else
            LWARNING("n=" << n << ":  Result GPU: " << resultCL << ", Result CPU: " << resultCPU << ", Relative Error: " << relError << " (FAILED)");
        numTests++;

        n = (1<<i) + 17;
        resultCL = voreenBlasCL.sNRM2(n, x);
        resultCPU = voreenBlasCPU.sNRM2(n, x);
        relError = std::abs((resultCL / resultCPU) - 1.f);
        if (relError < std::min(ERROR_THRESH*n, 0.1f)) {
            LINFO("n=" << n << ":  Result GPU: " << resultCL << ", Result CPU: " << resultCPU << ", Relative Error: " << relError);
            numSucc++;
        }
        else
            LWARNING("n=" << n << ":  Result GPU: " << resultCL << ", Result CPU: " << resultCPU << ", Relative Error: " << relError << " (FAILED)");
        numTests++;

    }

    delete[] x;

    if (numSucc == numTests) {
        LINFO("sNRM2 Randomized Tests (passed/overall): " << numSucc << "/" << numTests);
        return true;
    }
    else {
        LWARNING("sNRM2 Randomized Tests (passed/overall): " << numSucc << "/" << numTests);
        return false;
    }
}

bool testSSpMVEll() {

    float thresh = 5e-6;
    int numTests = 0;
    int numSucc = 0;

    LINFO("sSpMV-ELL Ground Truth Test: rows=8 cols=6 colsPerRow=3");

    EllpackMatrix<float> mat(8, 6, 3);
    mat.initializeBuffers();
    mat.setValue(0, 1, 2.f);
    mat.setValue(0, 3, -3.f);
    mat.setValue(0, 4, 5.f);
    mat.setValue(1, 0, -1.5);
    mat.setValue(1, 1, 9.f);
    mat.setValue(2, 1, 2.f);
    mat.setValue(2, 3, -3.5f);
    mat.setValue(2, 4, 7.f);
    mat.setValue(4, 0, 5.5f);
    mat.setValue(4, 5, -13.f);
    mat.setValue(5, 1, 8.f);
    mat.setValue(5, 2, 9.f);
    mat.setValue(5, 3, -5.f);
    mat.setValue(6, 0, 2.f);
    mat.setValue(6, 4, -3.f);
    mat.setValue(6, 5, -0.5f);
    mat.setValue(7, 0, 1.f);
    mat.setValue(7, 2, 2.f);
    mat.setValue(7, 5, 3.f);

    float* x = new float[6];
    x[0] = 1.f;
    x[1] = -3.f;
    x[2] = 0.f;
    x[3] = 5.f;
    x[4] = 0.5f;
    x[5] = 10.f;

    float* expRes = new float[8];
    expRes[0] = -18.5f;
    expRes[1] = -28.5;
    expRes[2] = -20.f;
    expRes[3] = 0.f;
    expRes[4] = -124.5f;
    expRes[5] = -49;
    expRes[6] = -4.5;
    expRes[7] = 31.f;

    float* yCL = new float[8];
    float* yCPU = new float[8];
    voreenBlasCL.sSpMVEll(mat, x, yCL);
    voreenBlasCPU.sSpMVEll(mat, x, yCPU);

    LINFO("Ground Truth: "  << expRes[0] << " " << expRes[1] << " " << expRes[2] << " " << expRes[3] << " " << expRes[4] << " " << expRes[5] << " " << expRes[6] << " " << expRes[7]);
    if (cmpVectors(8, expRes, yCPU, thresh)) {
        LINFO("Result CPU : " << yCPU[0]   << " " << yCPU[1]   << " " << yCPU[2]   << " " << yCPU[3]   << " " << yCPU[4]   << " " << yCPU[5]   << " " << yCPU[6]   << " " << yCPU[7] << " (passed)");
        numSucc++;
    }
    else
        LWARNING("Result CPU : " << yCPU[0]   << " " << yCPU[1]   << " " << yCPU[2]   << " " << yCPU[3]   << " " << yCPU[4]   << " " << yCPU[5]   << " " << yCPU[6]   << " " << yCPU[7] << " (FAILED)");
    if (cmpVectors(8, expRes, yCL, thresh)) {
        LINFO("Result OCL : " << yCL[0]    << " " << yCL[1]    << " " << yCL[2]    << " " << yCL[3]    << " " << yCL[4]    << " " << yCL[5]    << " " << yCL[6]    << " " << yCL[7] << " (passed)");
        numSucc++;
    }
    else
        LWARNING("Result OCL : " << yCL[0]    << " " << yCL[1]    << " " << yCL[2]    << " " << yCL[3]    << " " << yCL[4]    << " " << yCL[5]    << " " << yCL[6]    << " " << yCL[7] << " (FAILED)");
    numTests += 2;

    delete[] x;
    delete[] yCL;
    delete[] yCPU;
    delete[] expRes;


    LINFO("sSpMV-ELL Randomized Tests:");
    float* diffBuff = 0;
    for (int i=2; i<=MAX_MATRIX_SIZE; ++i) {

        {
            int numRows = (1<<i);
            int numCols = numRows;
            int numColsPerRow = 1;

            EllpackMatrix<float> randMat(numRows, numCols, numColsPerRow);
            randMat.initializeBuffers();
            randomizeEll(randMat);

            x = new float[randMat.getNumCols()];
            yCL = new float[randMat.getNumRows()];
            yCPU = new float[randMat.getNumRows()];
            diffBuff = new float[randMat.getNumRows()];
            for (size_t j=0; j<randMat.getNumCols(); ++j)
                x[j] = (static_cast<float>(rand()) / RAND_MAX) * 10.f - 5.f;

            voreenBlasCL.sSpMVEll(randMat, x, yCL);
            voreenBlasCPU.sSpMVEll(randMat, x, yCPU);

            voreenBlasCPU.sAXPY(numRows, yCL, yCPU, -1.f, diffBuff);
            float norm = voreenBlasCPU.sNRM2(numRows, diffBuff);
            if (norm < thresh*numRows) {
                LINFO("rows=" << randMat.getNumRows() << " cols=" << randMat.getNumCols() << " colsPerRow= " << randMat.getNumColsPerRow() << ": passed (" << norm << ")");
                numSucc++;
            }
            else {
                LWARNING("rows=" << randMat.getNumRows() << " cols=" << randMat.getNumCols() << " colsPerRow= " << randMat.getNumColsPerRow() << ": FAILED (" << norm << ")");
            }
            numTests++;

            delete[] x;
            delete[] yCL;
            delete[] yCPU;
            delete[] diffBuff;
        }

        {
            int numRows = (1<<i) - 1;
            int numCols = tgt::iround(numRows / 1.5) + 7;
            float maxCols;
            if (numRows < 1<<12)
                maxCols = 128.f;
            else if (numRows < 1<<16)
                maxCols = 32.f;
            else if (numRows < 1<<20)
                maxCols = 8.f;
            else
                maxCols = 4.f;
            int numColsPerRow = tgt::clamp(tgt::iround((static_cast<float>(rand()) / RAND_MAX)*maxCols), 2, numCols);

            EllpackMatrix<float> randMat(numRows, numCols, numColsPerRow);
            randMat.initializeBuffers();
            randomizeEll(randMat);

            x = new float[randMat.getNumCols()];
            yCL = new float[randMat.getNumRows()];
            yCPU = new float[randMat.getNumRows()];
            diffBuff = new float[randMat.getNumRows()];
            for (size_t j=0; j<randMat.getNumCols(); ++j)
                x[j] = (static_cast<float>(rand()) / RAND_MAX) * 10.f - 5.f;

            voreenBlasCL.sSpMVEll(randMat, x, yCL);
            voreenBlasCPU.sSpMVEll(randMat, x, yCPU);
            voreenBlasCPU.sAXPY(numRows, yCL, yCPU, -1.f, diffBuff);
            float norm = voreenBlasCPU.sNRM2(numRows, diffBuff);
            if (norm < thresh*numRows) {
                LINFO("rows=" << randMat.getNumRows() << " cols=" << randMat.getNumCols() << " colsPerRow= " << randMat.getNumColsPerRow() << ": passed (" << norm << ")");
                numSucc++;
            }
            else {
                LWARNING("rows=" << randMat.getNumRows() << " cols=" << randMat.getNumCols() << " colsPerRow= " << randMat.getNumColsPerRow() << ": FAILED (" << norm << ")");
            }
            numTests++;

            delete[] x;
            delete[] yCL;
            delete[] yCPU;
            delete[] diffBuff;
        }

        {
            int numRows = (1<<i) + 37;
            int numCols = tgt::iround(numRows*1.1f) + 3;
            float maxCols;
            if (numRows < 1<<12)
                maxCols = 128.f;
            else if (numRows < 1<<16)
                maxCols = 32.f;
            else if (numRows < 1<<20)
                maxCols = 8.f;
            else
                maxCols = 4.f;
            int numColsPerRow = tgt::clamp(tgt::iround((static_cast<float>(rand()) / RAND_MAX)*maxCols), 2, numCols);

            EllpackMatrix<float> randMat(numRows, numCols, numColsPerRow);
            randMat.initializeBuffers();
            randomizeEll(randMat);

            x = new float[randMat.getNumCols()];
            yCL = new float[randMat.getNumRows()];
            yCPU = new float[randMat.getNumRows()];
            diffBuff = new float[randMat.getNumRows()];
            for (size_t j=0; j<randMat.getNumCols(); ++j)
                x[j] = (static_cast<float>(rand()) / RAND_MAX) * 10.f - 5.f;

            voreenBlasCL.sSpMVEll(randMat, x, yCL);
            voreenBlasCPU.sSpMVEll(randMat, x, yCPU);
            voreenBlasCPU.sAXPY(numRows, yCL, yCPU, -1.f, diffBuff);
            float norm = voreenBlasCPU.sNRM2(numRows, diffBuff);
            if (norm < thresh*numRows) {
                LINFO("rows=" << randMat.getNumRows() << " cols=" << randMat.getNumCols() << " colsPerRow=" << randMat.getNumColsPerRow() << ": passed (" << norm << ")");
                numSucc++;
            }
            else {
                LWARNING("rows=" << randMat.getNumRows() << " cols=" << randMat.getNumCols() << " colsPerRow=" << randMat.getNumColsPerRow() << ": FAILED (" << norm << ")");
            }
            numTests++;

            delete[] x;
            delete[] yCL;
            delete[] yCPU;
            delete[] diffBuff;
        }
    }
    if (numSucc == numTests) {
        LINFO("sSpMV-ELL Randomized Tests (passed/overall): " << numSucc << "/" << numTests);
        return true;
    }
    else {
        LWARNING("sSpMV-ELL Randomized Tests (passed/overall): " << numSucc << "/" << numTests);
        return false;
    }
}

bool testHSpMVEll() {

    float thresh = 1e-5;
    int numTests = 0;
    int numSucc = 0;

    LINFO("hSpMV-ELL Ground Truth Test: rows=8 cols=6 colsPerRow=3");

    EllpackMatrix<int16_t> mat(8, 6, 3);
    mat.initializeBuffers();
    mat.setValue(0, 1, 4);
    mat.setValue(0, 3, -6);
    mat.setValue(0, 4, 10);
    mat.setValue(1, 0, -3);
    mat.setValue(1, 1, 18);
    mat.setValue(2, 1, 4);
    mat.setValue(2, 3, -7);
    mat.setValue(2, 4, 14);
    mat.setValue(4, 0, 11);
    mat.setValue(4, 5, -26);
    mat.setValue(5, 1, 16);
    mat.setValue(5, 2, 18);
    mat.setValue(5, 3, -10);
    mat.setValue(6, 0, 4);
    mat.setValue(6, 4, -6);
    mat.setValue(6, 5, -1);
    mat.setValue(7, 0, 2);
    mat.setValue(7, 2, 4);
    mat.setValue(7, 5, 6);

    float* x = new float[6];
    x[0] = 1.f;
    x[1] = -3.f;
    x[2] = 0.f;
    x[3] = 5.f;
    x[4] = 0.5f;
    x[5] = 10.f;

    float scale = 2.f / ((1<<15) - 1);
    float* expRes = new float[8];
    expRes[0] = -18.5f * scale;
    expRes[1] = -28.5f * scale;
    expRes[2] = -20.f * scale;
    expRes[3] = 0.f * scale;
    expRes[4] = -124.5f * scale;
    expRes[5] = -49.f * scale;
    expRes[6] = -4.5f * scale;
    expRes[7] = 31.f * scale;

    float* yCL = new float[8];
    float* yCPU = new float[8];
    voreenBlasCL.hSpMVEll(mat, x, yCL);
    voreenBlasCPU.hSpMVEll(mat, x, yCPU);

    LINFO("Ground Truth: "  << expRes[0] << " " << expRes[1] << " " << expRes[2] << " " << expRes[3] << " " << expRes[4] << " " << expRes[5] << " " << expRes[6] << " " << expRes[7]);
    if (cmpVectors(8, expRes, yCPU, thresh)) {
        LINFO("Result CPU : " << yCPU[0]   << " " << yCPU[1]   << " " << yCPU[2]   << " " << yCPU[3]   << " " << yCPU[4]   << " " << yCPU[5]   << " " << yCPU[6]   << " " << yCPU[7] << " (passed)");
        numSucc++;
    }
    else
        LWARNING("Result CPU : " << yCPU[0]   << " " << yCPU[1]   << " " << yCPU[2]   << " " << yCPU[3]   << " " << yCPU[4]   << " " << yCPU[5]   << " " << yCPU[6]   << " " << yCPU[7] << " (FAILED)");
    if (cmpVectors(8, expRes, yCL, thresh)) {
        LINFO("Result OCL : " << yCL[0]    << " " << yCL[1]    << " " << yCL[2]    << " " << yCL[3]    << " " << yCL[4]    << " " << yCL[5]    << " " << yCL[6]    << " " << yCL[7] << " (passed)");
        numSucc++;
    }
    else
        LWARNING("Result OCL : " << yCL[0]    << " " << yCL[1]    << " " << yCL[2]    << " " << yCL[3]    << " " << yCL[4]    << " " << yCL[5]    << " " << yCL[6]    << " " << yCL[7] << " (FAILED)");

    numTests += 2;

    delete[] x;
    delete[] yCL;
    delete[] yCPU;
    delete[] expRes;

    LINFO("hSpMV-ELL Randomized Tests:");
    float* diffBuff = 0;
    for (int i=2; i<=MAX_MATRIX_SIZE; ++i) {

        {
            int numRows = (1<<i);
            int numCols = numRows;
            int numColsPerRow = 1;

            EllpackMatrix<int16_t> randMat(numRows, numCols, numColsPerRow);
            randMat.initializeBuffers();
            randomizeEll(randMat);

            x = new float[randMat.getNumCols()];
            yCL = new float[randMat.getNumRows()];
            yCPU = new float[randMat.getNumRows()];
            diffBuff = new float[randMat.getNumRows()];
            for (size_t j=0; j<randMat.getNumCols(); ++j)
                x[j] = (static_cast<float>(rand()) / RAND_MAX) * 10.f - 5.f;

            voreenBlasCL.hSpMVEll(randMat, x, yCL);
            voreenBlasCPU.hSpMVEll(randMat, x, yCPU);

            voreenBlasCPU.sAXPY(numRows, yCL, yCPU, -1.f, diffBuff);
            float norm = voreenBlasCPU.sNRM2(numRows, diffBuff);
            if (norm < thresh*numRows) {
                LINFO("rows=" << randMat.getNumRows() << " cols=" << randMat.getNumCols() << " colsPerRow= " << randMat.getNumColsPerRow() << ": passed (" << norm << ")");
                numSucc++;
            }
            else {
                LWARNING("rows=" << randMat.getNumRows() << " cols=" << randMat.getNumCols() << " colsPerRow= " << randMat.getNumColsPerRow() << ": FAILED (" << norm << ")");
            }
            numTests++;

            delete[] x;
            delete[] yCL;
            delete[] yCPU;
            delete[] diffBuff;
        }

        {
            int numRows = (1<<i) - 1;
            int numCols = tgt::iround(numRows / 1.5) + 7;
            float maxCols;
            if (numRows < 1<<12)
                maxCols = 128.f;
            else if (numRows < 1<<16)
                maxCols = 32.f;
            else if (numRows < 1<<20)
                maxCols = 8.f;
            else
                maxCols = 4.f;
            int numColsPerRow = tgt::clamp(tgt::iround((static_cast<float>(rand()) / RAND_MAX)*maxCols), 2, numCols);

            EllpackMatrix<int16_t> randMat(numRows, numCols, numColsPerRow);
            randMat.initializeBuffers();
            randomizeEll(randMat);

            x = new float[randMat.getNumCols()];
            yCL = new float[randMat.getNumRows()];
            yCPU = new float[randMat.getNumRows()];
            diffBuff = new float[randMat.getNumRows()];
            for (size_t j=0; j<randMat.getNumCols(); ++j)
                x[j] = (static_cast<float>(rand()) / RAND_MAX) * 10.f - 5.f;

            voreenBlasCL.hSpMVEll(randMat, x, yCL);
            voreenBlasCPU.hSpMVEll(randMat, x, yCPU);
            voreenBlasCPU.sAXPY(numRows, yCL, yCPU, -1.f, diffBuff);
            float norm = voreenBlasCPU.sNRM2(numRows, diffBuff);
            if (norm < thresh*numRows) {
                LINFO("rows=" << randMat.getNumRows() << " cols=" << randMat.getNumCols() << " colsPerRow= " << randMat.getNumColsPerRow() << ": passed (" << norm << ")");
                numSucc++;
            }
            else {
                LWARNING("rows=" << randMat.getNumRows() << " cols=" << randMat.getNumCols() << " colsPerRow= " << randMat.getNumColsPerRow() << ": FAILED (" << norm << ")");
            }
            numTests++;

            delete[] x;
            delete[] yCL;
            delete[] yCPU;
            delete[] diffBuff;
        }

        {
            int numRows = (1<<i) + 37;
            int numCols = tgt::iround(numRows*1.1f) + 3;
            float maxCols;
            if (numRows < 1<<12)
                maxCols = 128.f;
            else if (numRows < 1<<16)
                maxCols = 32.f;
            else if (numRows < 1<<20)
                maxCols = 8.f;
            else
                maxCols = 4.f;
            int numColsPerRow = tgt::clamp(tgt::iround((static_cast<float>(rand()) / RAND_MAX)*maxCols), 2, numCols);

            EllpackMatrix<int16_t> randMat(numRows, numCols, numColsPerRow);
            randMat.initializeBuffers();
            randomizeEll(randMat);

            x = new float[randMat.getNumCols()];
            yCL = new float[randMat.getNumRows()];
            yCPU = new float[randMat.getNumRows()];
            diffBuff = new float[randMat.getNumRows()];
            for (size_t j=0; j<randMat.getNumCols(); ++j)
                x[j] = (static_cast<float>(rand()) / RAND_MAX) * 10.f - 5.f;

            voreenBlasCL.hSpMVEll(randMat, x, yCL);
            voreenBlasCPU.hSpMVEll(randMat, x, yCPU);
            voreenBlasCPU.sAXPY(numRows, yCL, yCPU, -1.f, diffBuff);
            float norm = voreenBlasCPU.sNRM2(numRows, diffBuff);
            if (norm < thresh*numRows) {
                LINFO("rows=" << randMat.getNumRows() << " cols=" << randMat.getNumCols() << " colsPerRow=" << randMat.getNumColsPerRow() << ": passed (" << norm << ")");
                numSucc++;
            }
            else {
                LWARNING("rows=" << randMat.getNumRows() << " cols=" << randMat.getNumCols() << " colsPerRow=" << randMat.getNumColsPerRow() << ": FAILED (" << norm << ")");
            }
            numTests++;

            delete[] x;
            delete[] yCL;
            delete[] yCPU;
            delete[] diffBuff;
        }
    }
    if (numSucc == numTests) {
        LINFO("hSpMV-ELL Randomized Tests (passed/overall): " << numSucc << "/" << numTests);
        return true;
    }
    else {
        LWARNING("hSpMV-ELL Randomized Tests (passed/overall): " << numSucc << "/" << numTests);
        return false;
    }
}

bool testSSpInnerProductEll() {

    LINFO("sInnerProduct-ELL Randomized Tests:");
    int numTests = 0;
    int numSucc = 0;
    for (int i=5; i<=MAX_MATRIX_SIZE; ++i) {

        {
            numTests++;

            int numRows = (1<<i);
            int numColsPerRow = std::min(7, numRows);

            EllpackMatrix<float> randMat(numRows, numRows, numColsPerRow);
            randMat.initializeBuffers();
            randomizeEllPositiveDefinite(randMat);
            if (!randMat.isSymmetric()) {
                LWARNING("rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": FAILED (matrix not symmetric)");
                continue;
            }

            float* x = new float[numRows];
            float* y = new float[numRows];
            for (int j=0; j<numRows; ++j) {
                x[j] = (static_cast<float>(rand()) / RAND_MAX) * 2.f - 0.5f;
                y[j] = (static_cast<float>(rand()) / RAND_MAX) * 2.f - 0.5f;
            }

            float resultCL = voreenBlasCL.sSpInnerProductEll(randMat, x, y);
            float resultCPU = voreenBlasCPU.sSpInnerProductEll(randMat, x, y);

            float relError = std::abs((resultCL / resultCPU) - 1.f);
            if (relError < std::min(ERROR_THRESH*numRows, 0.1f)) {
                LINFO("(pos. definite) n=" << numRows << ":  Result GPU: " << resultCL << ", Result CPU: " << resultCPU << ", Relative Error: " << relError);
                numSucc++;
            }
            else
                LWARNING("(pos. definite) n=" << numRows << ":  Result GPU: " << resultCL << ", Result CPU: " << resultCPU << ", Relative Error: " << relError << " (FAILED)");

            delete[] x;
            delete[] y;
        }

        {
            numTests++;

            int numRows = (1<<i) - 1;
            int numColsPerRow = std::min(7, numRows);

            EllpackMatrix<float> randMat(numRows, numRows, numColsPerRow);
            randMat.initializeBuffers();
            randomizeEll(randMat);

            float* x = new float[numRows];
            float* y = new float[numRows];
            for (int j=0; j<numRows; ++j) {
                x[j] = (static_cast<float>(rand()) / RAND_MAX) * 2.f - 0.5f;
                y[j] = (static_cast<float>(rand()) / RAND_MAX) * 2.f - 0.5f;
            }

            float resultCL = voreenBlasCL.sSpInnerProductEll(randMat, x, y);
            float resultCPU = voreenBlasCPU.sSpInnerProductEll(randMat, x, y);

            float relError = std::abs((resultCL / resultCPU) - 1.f);
            if (relError < std::min(ERROR_THRESH*numRows, 0.1f)) {
                LINFO("n=" << numRows << ":  Result GPU: " << resultCL << ", Result CPU: " << resultCPU << ", Relative Error: " << relError);
                numSucc++;
            }
            else
                LWARNING("n=" << numRows << ":  Result GPU: " << resultCL << ", Result CPU: " << resultCPU << ", Relative Error: " << relError << " (FAILED)");

            delete[] x;
            delete[] y;
        }

        {
            numTests++;

            int numRows = (1<<i) + 17;
            int numColsPerRow = std::min(7, numRows);

            EllpackMatrix<float> randMat(numRows, numRows, numColsPerRow);
            randMat.initializeBuffers();
            randomizeEll(randMat);

            float* x = new float[numRows];
            float* y = new float[numRows];
            for (int j=0; j<numRows; ++j) {
                x[j] = (static_cast<float>(rand()) / RAND_MAX) * 2.f - 0.5f;
                y[j] = (static_cast<float>(rand()) / RAND_MAX) * 2.f - 0.5f;
            }

            float resultCL = voreenBlasCL.sSpInnerProductEll(randMat, x, y);
            float resultCPU = voreenBlasCPU.sSpInnerProductEll(randMat, x, y);

            float relError = std::abs((resultCL / resultCPU) - 1.f);
            if (relError < std::min(ERROR_THRESH*numRows, 0.1f)) {
                LINFO("n=" << numRows << ":  Result GPU: " << resultCL << ", Result CPU: " << resultCPU << ", Relative Error: " << relError);
                numSucc++;
            }
            else
                LWARNING("n=" << numRows << ":  Result GPU: " << resultCL << ", Result CPU: " << resultCPU << ", Relative Error: " << relError << " (FAILED)");

            delete[] x;
            delete[] y;
        }

    }

    if (numSucc == numTests) {
        LINFO("sInnerProduct-ELL Randomized Tests (passed/overall): " << numSucc << "/" << numTests);
        return true;
    }
    else {
        LWARNING("sInnerProduct-ELL Randomized Tests Tests (passed/overall): " << numSucc << "/" << numTests);
        return false;
    }
}
struct NoopProgressReporter : public ProgressReporter {
    virtual void setProgressMessage(const std::string& message) { }
    virtual std::string getProgressMessage() const {return std::string();}
};

bool testSSpConjGradEll() {

    NoopProgressReporter progress;

    LINFO("sSpConjGrad-ELL Randomized Tests:");
    int numTests = 0;
    int numSucc = 0;
    for (int i=3; i<=MAX_MATRIX_SIZE; ++i) {

        {
            int numRows = (1<<i);
            int numColsPerRow = std::min(7, numRows);

            // initialize test matrix
            EllpackMatrix<float> randMat(numRows, numRows, numColsPerRow);
            randMat.initializeBuffers();
            randomizeEllPositiveDefinite(randMat);
            if (!randMat.isSymmetric()) {
                LWARNING("rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": FAILED (matrix not symmetric)");
                continue;
            }

            // initialize vectors
            float* x = new float[numRows];
            float* y = new float[numRows];
            float* tempBuf = new float[numRows];
            for (int j=0; j<numRows; ++j)
                x[j] = (static_cast<float>(rand()) / RAND_MAX) * 1.f;

            // compute GPU solution
            clock_t start, end;
            start = clock();
            int iterations = voreenBlasCL.sSpConjGradEll(randMat, x, y, 0, VoreenBlas::NoPreconditioner, numRows*1e-8f, 1000, progress);
            end = clock();
            float duration = static_cast<float>(end-start)/CLOCKS_PER_SEC;

            // check result
            voreenBlasCPU.sSpMVEll(randMat, y, tempBuf);
            voreenBlasCPU.sAXPY(numRows, x, tempBuf, -1.f, tempBuf);
            float normGPU = voreenBlasCPU.sNRM2(numRows, tempBuf);
            if (normGPU < numRows*ERROR_THRESH) {
                LINFO("GPU: rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": passed (error: " << normGPU << ", iterations: " << iterations << ", duration: " << duration << " sec)");
                numSucc++;
            }
            else {
                LWARNING("GPU: rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": FAILED (" << normGPU << ")");
            }
            numTests++;

            // CPU solution
            if (i <= 18) {
                start = clock();
                iterations = voreenBlasCPU.sSpConjGradEll(randMat, x, y, 0, VoreenBlas::NoPreconditioner, numRows*1e-8f, 1000, progress);
                end = clock();
                duration = static_cast<float>(end-start)/CLOCKS_PER_SEC;

                voreenBlasCPU.sSpMVEll(randMat, y, tempBuf);
                voreenBlasCPU.sAXPY(numRows, x, tempBuf, -1.f, tempBuf);
                float normCPU = voreenBlasCPU.sNRM2(numRows, tempBuf);
                if (normCPU < numRows*ERROR_THRESH) {
                    LINFO("CPU: rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": passed (error: " << normCPU << ", iterations: " << iterations << ", duration: " << duration << " sec)");
                    numSucc++;
                }
                else {
                    LWARNING("CPU: rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": FAILED (" << normCPU << ")");
                }
                numTests++;
            }

            delete[] x;
            delete[] y;
            delete[] tempBuf;
        }

        {
            int numRows = tgt::iround((1<<i) * 0.7f);
            int numColsPerRow = std::min(6, numRows);

            // initialize test matrix
            EllpackMatrix<float> randMat(numRows, numRows, numColsPerRow);
            randMat.initializeBuffers();
            randomizeEllPositiveDefinite(randMat);
            if (!randMat.isSymmetric()) {
                LWARNING("rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": FAILED (matrix not symmetric)");
                continue;
            }

            // initialize vectors
            float* x = new float[numRows];
            float* y = new float[numRows];
            float* tempBuf = new float[numRows];
            for (int j=0; j<numRows; ++j)
                x[j] = (static_cast<float>(rand()) / RAND_MAX) * 2.f - 0.5f;

            // compute GPU solution
            clock_t start, end;
            start = clock();
            int iterations = voreenBlasCL.sSpConjGradEll(randMat, x, y, 0, VoreenBlas::NoPreconditioner, numRows*1e-8f, 1000, progress);
            end = clock();
            float duration = static_cast<float>(end-start)/CLOCKS_PER_SEC;

            // check result
            voreenBlasCPU.sSpMVEll(randMat, y, tempBuf);
            voreenBlasCPU.sAXPY(numRows, x, tempBuf, -1.f, tempBuf);
            float normGPU = voreenBlasCPU.sNRM2(numRows, tempBuf);
            if (normGPU < numRows*ERROR_THRESH) {
                LINFO("GPU: rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": passed (error: " << normGPU << ", iterations: " << iterations << ", duration: " << duration << " sec)");
                numSucc++;
            }
            else {
                LWARNING("GPU: rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": FAILED (" << normGPU << ")");
            }
            numTests++;

            // CPU solution
            if (i <= 18) {
                start = clock();
                iterations = voreenBlasCPU.sSpConjGradEll(randMat, x, y, 0, VoreenBlas::NoPreconditioner, numRows*1e-8f, 1000, progress);
                end = clock();
                duration = static_cast<float>(end-start)/CLOCKS_PER_SEC;

                voreenBlasCPU.sSpMVEll(randMat, y, tempBuf);
                voreenBlasCPU.sAXPY(numRows, x, tempBuf, -1.f, tempBuf);
                float normCPU = voreenBlasCPU.sNRM2(numRows, tempBuf);
                if (normCPU < numRows*ERROR_THRESH) {
                    LINFO("CPU: rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": passed (error: " << normCPU << ", iterations: " << iterations << ", duration: " << duration << " sec)");
                    numSucc++;
                }
                else {
                    LWARNING("CPU: rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": FAILED (" << normCPU << ")");
                }
                numTests++;
            }

            delete[] x;
            delete[] y;
            delete[] tempBuf;
        }

    }

    if (numSucc == numTests) {
        LINFO("sSpConjGrad-ELL Randomized Tests (passed/overall): " << numSucc << "/" << numTests);
        return true;
    }
    else {
        LWARNING("sSpConjGrad-ELL Randomized Tests (passed/overall): " << numSucc << "/" << numTests);
        return false;
    }
}


bool testHSpConjGradEll() {

    LINFO("hSpConjGrad-ELL Randomized Tests:");
    int numTests = 0;
    int numSucc = 0;
    for (int i=3; i<=MAX_MATRIX_SIZE; ++i) {

        {
            int numRows = (1<<i);
            int numColsPerRow = std::min(7, numRows);

            // initialize test matrix
            EllpackMatrix<int16_t> randMat(numRows, numRows, numColsPerRow);
            randMat.initializeBuffers();
            randomizeEllPositiveDefinite(randMat);
            if (!randMat.isSymmetric()) {
                LWARNING("rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": FAILED (matrix not symmetric)");
                continue;
            }

            // initialize vectors
            float* x = new float[numRows];
            float* y = new float[numRows];
            float* tempBuf = new float[numRows];
            for (int j=0; j<numRows; ++j)
                x[j] = (static_cast<float>(rand()) / RAND_MAX) * 1.f;

            // compute GPU solution
            clock_t start, end;
            start = clock();
            int iterations = voreenBlasCL.hSpConjGradEll(randMat, x, y, 0, numRows*1e-8f);
            end = clock();
            float duration = static_cast<float>(end-start)/CLOCKS_PER_SEC;

            // check result
            voreenBlasCPU.hSpMVEll(randMat, y, tempBuf);
            voreenBlasCPU.sAXPY(numRows, x, tempBuf, -1.f, tempBuf);
            float normGPU = voreenBlasCPU.sNRM2(numRows, tempBuf);
            if (normGPU < numRows*ERROR_THRESH) {
                LINFO("GPU: rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": passed (error: " << normGPU << ", iterations: " << iterations << ", duration: " << duration << " sec)");
                numSucc++;
            }
            else {
                LWARNING("GPU: rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": FAILED (" << normGPU << ")");
            }
            numTests++;

            // CPU solution
            if (i <= 18) {
                start = clock();
                iterations = voreenBlasCPU.hSpConjGradEll(randMat, x, y, 0, numRows*1e-8f);
                end = clock();
                duration = static_cast<float>(end-start)/CLOCKS_PER_SEC;

                voreenBlasCPU.hSpMVEll(randMat, y, tempBuf);
                voreenBlasCPU.sAXPY(numRows, x, tempBuf, -1.f, tempBuf);
                float normCPU = voreenBlasCPU.sNRM2(numRows, tempBuf);
                if (normCPU < numRows*ERROR_THRESH) {
                    LINFO("CPU: rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": passed (error: " << normCPU << ", iterations: " << iterations << ", duration: " << duration << " sec)");
                    numSucc++;
                }
                else {
                    LWARNING("CPU: rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": FAILED (" << normCPU << ")");
                }
                numTests++;
            }

            delete[] x;
            delete[] y;
            delete[] tempBuf;
        }

        {
            int numRows = tgt::iround((1<<i) * 0.7f);
            int numColsPerRow = std::min(6, numRows);

            // initialize test matrix
            EllpackMatrix<int16_t> randMat(numRows, numRows, numColsPerRow);
            randMat.initializeBuffers();
            randomizeEllPositiveDefinite(randMat);
            if (!randMat.isSymmetric()) {
                LWARNING("rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": FAILED (matrix not symmetric)");
                continue;
            }

            // initialize vectors
            float* x = new float[numRows];
            float* y = new float[numRows];
            float* tempBuf = new float[numRows];
            for (int j=0; j<numRows; ++j)
                x[j] = (static_cast<float>(rand()) / RAND_MAX) * 2.f - 0.5f;

            // compute GPU solution
            clock_t start, end;
            start = clock();
            int iterations = voreenBlasCL.hSpConjGradEll(randMat, x, y, 0, numRows*1e-8f);
            end = clock();
            float duration = static_cast<float>(end-start)/CLOCKS_PER_SEC;

            // check result
            voreenBlasCPU.hSpMVEll(randMat, y, tempBuf);
            voreenBlasCPU.sAXPY(numRows, x, tempBuf, -1.f, tempBuf);
            float normGPU = voreenBlasCPU.sNRM2(numRows, tempBuf);
            if (normGPU < numRows*ERROR_THRESH) {
                LINFO("GPU: rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": passed (error: " << normGPU << ", iterations: " << iterations << ", duration: " << duration << " sec)");
                numSucc++;
            }
            else {
                LWARNING("GPU: rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": FAILED (" << normGPU << ")");
            }
            numTests++;

            // CPU solution
            if (i <= 18) {
                start = clock();
                iterations = voreenBlasCPU.hSpConjGradEll(randMat, x, y, 0, numRows*1e-8f);
                end = clock();
                duration = static_cast<float>(end-start)/CLOCKS_PER_SEC;

                voreenBlasCPU.hSpMVEll(randMat, y, tempBuf);
                voreenBlasCPU.sAXPY(numRows, x, tempBuf, -1.f, tempBuf);
                float normCPU = voreenBlasCPU.sNRM2(numRows, tempBuf);
                if (normCPU < numRows*ERROR_THRESH) {
                    LINFO("CPU: rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": passed (error: " << normCPU << ", iterations: " << iterations << ", duration: " << duration << " sec)");
                    numSucc++;
                }
                else {
                    LWARNING("CPU: rows=" << randMat.getNumRows() << " colsPerRow= " << randMat.getNumColsPerRow() << ": FAILED (" << normCPU << ")");
                }
                numTests++;
            }

            delete[] x;
            delete[] y;
            delete[] tempBuf;
        }

    }

    if (numSucc == numTests) {
        LINFO("hSpConjGrad-ELL Randomized Tests (passed/overall): " << numSucc << "/" << numTests);
        return true;
    }
    else {
        LWARNING("hSpConjGrad-ELL Randomized Tests (passed/overall): " << numSucc << "/" << numTests);
        return false;
    }

}
void randomizeEll(EllpackMatrix<float>& mat) {

    tgt::vec2 range(-5.f, 10.f);

    for (size_t row=0; row < mat.getNumRows(); ++row) {
        // leave 1/100 rows empty
        if ((rand() % 100) == 0)
            continue;

        // fill 1/10 of the cols only partially
        size_t significantCols = tgt::iround(10.f * (static_cast<float>(rand()) / RAND_MAX) * mat.getNumColsPerRow() );
        significantCols = std::min(significantCols, mat.getNumColsPerRow());
        int lastCol = -1;
        for (size_t colIndex = 0; colIndex < significantCols; ++colIndex) {
            size_t col = (lastCol + 1) + tgt::iround( (static_cast<float>(rand()) / RAND_MAX) * (mat.getNumColsPerRow()-1 - lastCol));
            if (col < mat.getNumColsPerRow()) {
                float val = range.x + (static_cast<float>(rand()) / RAND_MAX) * (range.y - range.x);
                mat.setValueByIndex(row, col, colIndex, val);
                lastCol = (int)col;
            }
        }
    }
}

void randomizeEll(EllpackMatrix<int16_t>& mat) {

    tgt::Vector2<int16_t> range(-(1<<15) + 1, (1<<15) - 1);

    for (size_t row=0; row < mat.getNumRows(); ++row) {
        // leave 1/100 rows empty
        if ((rand() % 100) == 0)
            continue;

        // fill 1/10 of the cols only partially
        size_t significantCols = tgt::iround(10.f * (static_cast<float>(rand()) / RAND_MAX) * mat.getNumColsPerRow() );
        significantCols = std::min(significantCols, mat.getNumColsPerRow());
        int lastCol = -1;
        for (size_t colIndex = 0; colIndex < significantCols; ++colIndex) {
            size_t col = (lastCol + 1) + tgt::iround( (static_cast<float>(rand()) / RAND_MAX) * (mat.getNumColsPerRow()-1 - lastCol));
            if (col < mat.getNumColsPerRow()) {
                int16_t val = static_cast<int16_t>(range.x + (static_cast<float>(rand()) / RAND_MAX) * (range.y - range.x));
                mat.setValueByIndex(row, col, colIndex, val);
                lastCol = (int)col;
            }
        }
    }
}

void randomizeEllPositiveDefinite(EllpackMatrix<float>& mat) {

    tgtAssert(mat.isQuadratic(), "quadratic matrix expected");

    // generate random volume
    int volDim = tgt::iceil(powf((float)mat.getNumRows(), 1.f/3.f));
    tgt::ivec3 dim(volDim, volDim, tgt::iceil(volDim*1.1f));
    const int numVoxels =  dim.x*dim.y*dim.z;
    const int numSeeds = numVoxels - (int)mat.getNumRows();
    tgtAssert(numSeeds > 0, "Invalid seed count");

    uint8_t* volume = new uint8_t[numVoxels];
    bool* seeds = new bool[numVoxels];
    for (int i=0; i<numVoxels; i++) {
        volume[i] = tgt::ifloor((static_cast<float>(rand()) / RAND_MAX) * 255.f);
        seeds[i] = false;
    }

    int seedsGenerated = 0;
    while (seedsGenerated < numSeeds) {
        for (int i=0; i<numVoxels; i++) {
            if (seedsGenerated == numSeeds)
                break;
            if (!seeds[i] && ((rand() % 10) == 0)) {
                seeds[i] = true;
                seedsGenerated++;
            }
        }
    }

    int* volIndexToRow = new int[numVoxels];
    int curRow = 0;
    for (int i=0; i<numVoxels; i++) {
        if (!seeds[i]) {
            volIndexToRow[i] = curRow;
            curRow++;
        }
        else {
            volIndexToRow[i] = -1;
        }
    }


    float maxWeightSum = 0.f;
    tgt::ivec3 coords;
    for (int z=0; z<dim.z; z++) {
        for (int y=0; y<dim.y; y++) {
            for (int x=0; x<dim.x; x++) {

                int index = volumeCoordsToIndex(x, y, z, dim);
                if (seeds[index])
                    continue;

                size_t row = volIndexToRow[index];
                tgtAssert(row < mat.getNumCols(), "Invalid row");

                float weightSum = 0.f;
                int curIntensity = volume[index];

                // x-neighbors
                if (mat.getNumColsPerRow() >= 3) {
                    if (x > 0) {
                        int x0 = volumeCoordsToIndex(x-1, y, z, dim);
                        float weight = (std::abs(curIntensity - volume[x0]) + 1) / 255.f;
                        weightSum += weight;
                        if (!seeds[x0]) {
                            size_t rowX0 = volIndexToRow[x0];
                            tgtAssert(rowX0 < mat.getNumCols(), "Invalid row");
                            mat.setValue(row, rowX0, -weight);
                            mat.setValue(rowX0, row, -weight);
                        }
                    }
                    if (x < dim.x-1) {
                        size_t x1 = volumeCoordsToIndex(x+1, y, z, dim);
                        float weight = (std::abs(curIntensity - volume[x1]) + 1) / 255.f;
                        weightSum += weight;
                        if (!seeds[x1]) {
                            size_t rowX1 = volIndexToRow[x1];
                            tgtAssert(rowX1 < mat.getNumCols(), "Invalid row");
                            mat.setValue(row, rowX1, weight);
                            mat.setValue(rowX1, row, weight);
                        }
                    }
                }

                // y-neighbors
                if (mat.getNumColsPerRow() >= 5) {
                    if (y > 0) {
                        size_t y0 = volumeCoordsToIndex(x, y-1, z, dim);
                        float weight = std::abs(curIntensity - volume[y0]) / 255.f;
                        weightSum += weight;
                        if (!seeds[y0]) {
                            size_t rowY0 = volIndexToRow[y0];
                            tgtAssert(rowY0 < mat.getNumCols(), "Invalid row");
                            mat.setValue(row, rowY0, -weight);
                            mat.setValue(rowY0, row, -weight);
                        }
                    }
                    if (y < dim.y-1) {
                        int y1 = volumeCoordsToIndex(x, y+1, z, dim);
                        float weight = std::abs(curIntensity - volume[y1]) / 255.f;
                        weightSum += weight;
                        if (!seeds[y1]) {
                            size_t rowY1 = volIndexToRow[y1];
                            tgtAssert(rowY1 < mat.getNumCols(), "Invalid row");
                            mat.setValue(row, rowY1, weight);
                            mat.setValue(rowY1, row, weight);
                        }
                    }
                }

                // z-neighbors
                if (mat.getNumColsPerRow() >= 7) {
                    if (z > 0) {
                        size_t z0 = volumeCoordsToIndex(x, y, z-1, dim);
                        float weight = std::abs(curIntensity - volume[z0]) / 255.f;
                        weightSum += weight;
                        if (!seeds[z0]) {
                            size_t rowZ0 = volIndexToRow[z0];
                            tgtAssert(rowZ0 < mat.getNumCols(), "Invalid row");
                            mat.setValue(row, rowZ0, -weight);
                            mat.setValue(rowZ0, row, -weight);
                        }
                    }
                    if (z < dim.z-1) {
                        size_t z1 = volumeCoordsToIndex(x, y, z+1, dim);
                        float weight = std::abs(curIntensity - volume[z1]) / 255.f;
                        weightSum += weight;
                        if (!seeds[z1]) {
                            size_t rowZ1 = volIndexToRow[z1];
                            tgtAssert(rowZ1 < mat.getNumCols(), "Invalid row");
                            mat.setValue(row, rowZ1, weight);
                            mat.setValue(rowZ1, row, weight);
                        }
                    }
                }

                mat.setValue(row, row, weightSum);
                maxWeightSum = std::max(maxWeightSum, weightSum);
            }
        }
    }

    delete[] seeds;
    delete[] volume;
    delete[] volIndexToRow;
}

void randomizeEllPositiveDefinite(EllpackMatrix<int16_t>& mat) {
    tgtAssert(mat.isQuadratic(), "quadratic matrix expected");

    // generate random volume
    int volDim = tgt::iceil(powf((float)mat.getNumRows(), 1.f/3.f));
    tgt::ivec3 dim(volDim, volDim, tgt::iceil(volDim*1.1f));
    const int numVoxels =  dim.x*dim.y*dim.z;
    const int numSeeds = numVoxels - (int)mat.getNumRows();
    tgtAssert(numSeeds > 0, "Invalid seed count");

    uint8_t* volume = new uint8_t[numVoxels];
    bool* seeds = new bool[numVoxels];
    for (int i=0; i<numVoxels; i++) {
        volume[i] = tgt::ifloor((static_cast<float>(rand()) / RAND_MAX) * 255.f);
        seeds[i] = false;
    }

    int seedsGenerated = 0;
    while (seedsGenerated < numSeeds) {
        for (int i=0; i<numVoxels; i++) {
            if (seedsGenerated == numSeeds)
                break;
            if (!seeds[i] && ((rand() % 10) == 0)) {
                seeds[i] = true;
                seedsGenerated++;
            }
        }
    }

    int* volIndexToRow = new int[numVoxels];
    int curRow = 0;
    for (int i=0; i<numVoxels; i++) {
        if (!seeds[i]) {
            volIndexToRow[i] = curRow;
            curRow++;
        }
        else {
            volIndexToRow[i] = -1;
        }
    }


    int* diagonal = new int[mat.getNumCols()];
    float maxValue = static_cast<float>((1<<15) - 1);

    int maxWeightSum = 0;
    tgt::ivec3 coords;
    for (int z=0; z<dim.z; z++) {
        for (int y=0; y<dim.y; y++) {
            for (int x=0; x<dim.x; x++) {

                int index = volumeCoordsToIndex(x, y, z, dim);
                if (seeds[index])
                    continue;

                size_t row = volIndexToRow[index];
                tgtAssert(row < mat.getNumCols(), "Invalid row");

                int weightSum = 0;
                int16_t curIntensity = volume[index];

                // x-neighbors
                if (mat.getNumColsPerRow() >= 3) {
                    if (x > 0) {
                        size_t x0 = volumeCoordsToIndex(x-1, y, z, dim);
                        int16_t weight = std::abs(curIntensity - volume[x0]) + 1;
                        weightSum += (int)weight;
                        if (!seeds[x0]) {
                            size_t rowX0 = volIndexToRow[x0];
                            tgtAssert(rowX0 < mat.getNumCols(), "Invalid row");
                            mat.setValue(row, rowX0, -weight);
                            mat.setValue(rowX0, row, -weight);
                        }
                    }
                    if (x < dim.x-1) {
                        size_t x1 = volumeCoordsToIndex(x+1, y, z, dim);
                        int16_t weight = std::abs(curIntensity - volume[x1]) + 1;
                        weightSum += (int)weight;
                        if (!seeds[x1]) {
                            size_t rowX1 = volIndexToRow[x1];
                            tgtAssert(rowX1 < mat.getNumCols(), "Invalid row");
                            mat.setValue(row, rowX1, weight);
                            mat.setValue(rowX1, row, weight);
                        }
                    }
                }

                // y-neighbors
                if (mat.getNumColsPerRow() >= 5) {
                    if (y > 0) {
                        size_t y0 = volumeCoordsToIndex(x, y-1, z, dim);
                        int16_t weight = std::abs(curIntensity - volume[y0]);
                        weightSum += (int)weight;
                        if (!seeds[y0]) {
                            size_t rowY0 = volIndexToRow[y0];
                            tgtAssert(rowY0 < mat.getNumCols(), "Invalid row");
                            mat.setValue(row, rowY0, -weight);
                            mat.setValue(rowY0, row, -weight);
                        }
                    }
                    if (y < dim.y-1) {
                        size_t y1 = volumeCoordsToIndex(x, y+1, z, dim);
                        int16_t weight = std::abs(curIntensity - volume[y1]);
                        weightSum += (int)weight;
                        if (!seeds[y1]) {
                            size_t rowY1 = volIndexToRow[y1];
                            tgtAssert(rowY1 < mat.getNumCols(), "Invalid row");
                            mat.setValue(row, rowY1, weight);
                            mat.setValue(rowY1, row, weight);
                        }
                    }
                }

                // z-neighbors
                if (mat.getNumColsPerRow() >= 7) {
                    if (z > 0) {
                        size_t z0 = volumeCoordsToIndex(x, y, z-1, dim);
                        int16_t weight = std::abs(curIntensity - volume[z0]);
                        weightSum += (int)weight;
                        if (!seeds[z0]) {
                            size_t rowZ0 = volIndexToRow[z0];
                            tgtAssert(rowZ0 < mat.getNumCols(), "Invalid row");
                            mat.setValue(row, rowZ0, -weight);
                            mat.setValue(rowZ0, row, -weight);
                        }
                    }
                    if (z < dim.z-1) {
                        size_t z1 = volumeCoordsToIndex(x, y, z+1, dim);
                        int16_t weight = std::abs(curIntensity - volume[z1]);
                        weightSum += (int)weight;
                        if (!seeds[z1]) {
                            size_t rowZ1 = volIndexToRow[z1];
                            tgtAssert(rowZ1 < mat.getNumCols(), "Invalid row");
                            mat.setValue(row, rowZ1, weight);
                            mat.setValue(rowZ1, row, weight);
                        }
                    }
                }

                //mat.setValue(row, row, weightSum);
                diagonal[row] = weightSum;
                maxWeightSum = std::max(weightSum, maxWeightSum);
            }
        }
    }

    float scale = maxValue / maxWeightSum;
    for (size_t i=0; i<mat.getNumRows()*mat.getNumColsPerRow(); i++)
        mat.getMatrix()[i] = static_cast<int16_t>(mat.getMatrix()[i] * scale);

    for (size_t row=0; row<mat.getNumRows(); row++) {
        mat.setValue(row, row, static_cast<int16_t>(diagonal[row] * scale));
    }


    delete[] diagonal;
    delete[] seeds;
    delete[] volume;
    delete[] volIndexToRow;
}

bool cmpVectors(size_t vecSize, float* vecx, float* vecy, float relThresh) {
    for (size_t i=0; i<vecSize; ++i) {
        if ( ((std::abs(vecx[i]) > 1e2f*relThresh) || (std::abs(vecy[i]) > 1e2f*relThresh)) &&
             (std::abs((vecx[i] / vecy[i]) - 1.f) > relThresh) ) {
            return false;
        }
    }
    return true;
}
