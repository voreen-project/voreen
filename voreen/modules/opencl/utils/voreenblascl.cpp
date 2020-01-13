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

#include "voreenblascl.h"
#include "modules/opencl/openclmodule.h"

#include "voreen/core/voreenapplication.h"

namespace voreen {

using namespace cl;

const std::string VoreenBlasCL::loggerCat_("voreen.opencl.VoreenBlasCL");

VoreenBlasCL::VoreenBlasCL() :
    initialized_(false),
    opencl_(0),
    context_(0),
    //device_(0),
    queue_(0),
    prog_(0)
{
}

VoreenBlasCL::~VoreenBlasCL() {
    if (initialized_)
        deinitialize();
}

void VoreenBlasCL::initialize() {
    tgtAssert(!initialized_, "VoreenBlasCL has already been initialized");

    // acquire OpenCL resources
    if (!OpenCLModule::getInstance()->getCLContext())
        throw VoreenException("No OpenCL context created");

    opencl_ = OpenCLModule::getInstance()->getOpenCL();
    context_ = OpenCLModule::getInstance()->getCLContext();
    queue_ = OpenCLModule::getInstance()->getCLCommandQueue();
    device_ = OpenCLModule::getInstance()->getCLDevice();

    if (!opencl_ || !context_ || !queue_/* || !device_*/) {
        throw VoreenException("Failed to acquire OpenCL resources");
    }

    // load voreenblas.cl
    std::string kernelFile = VoreenApplication::app()->getModulePath("opencl") + "/cl/voreenblas.cl";
    LINFO("Loading program " << kernelFile);
    prog_ = OpenCLModule::getInstance()->loadProgram(kernelFile);

    initialized_ = true;
}

void VoreenBlasCL::deinitialize() {
    tgtAssert(initialized_, "VoreenBlasCL has not been intialized");

    delete prog_;
    prog_ = 0;

    // borrowed references
    queue_ = 0;
    context_ = 0;
    //device_ = 0;
    opencl_ = 0;

    initialized_ = false;
}

bool VoreenBlasCL::isInitialized() const {
    return initialized_;
}

void VoreenBlasCL::sAXPY(size_t vecSize, const float* vecx, const float* vecy, float alpha, float* result) const {

    if (!initialized_) {
        LERROR("Not initialized. Aborting.");
        return;
    }

    Kernel* kernel = prog_->getKernel("sAXPY");
    if (!kernel) {
        LERROR("No kernel 'sAXPY' found");
        return;
    }

    int workSize = 1024 << 2;

    Buffer xBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*vecSize, const_cast<float*>(vecx));
    Buffer yBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*vecSize, const_cast<float*>(vecy));

    Buffer resultBuffer(context_, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*vecSize);

    kernel->setArg(0, vecSize);
    kernel->setArg(1, xBuffer);
    kernel->setArg(2, yBuffer);
    kernel->setArg(3, alpha);
    kernel->setArg(4, resultBuffer);

    queue_->enqueue(kernel, workSize);

    queue_->enqueueReadBuffer(&resultBuffer, (void*)(result), true);
    queue_->finish();
}

float VoreenBlasCL::sDOT(size_t vecSize, const float* vecx, const float* vecy) const {

    if (!initialized_) {
        LERROR("Not initialized. Aborting.");
        return -1.f;
    }

    Kernel* kernel = prog_->getKernel("sDOT");
    if (!kernel) {
        LERROR("No kernel 'sDOT' found");
        return -1.f;
    }

    int workSize = 1024 << 2;
    int localWorkSize = std::min<int>(512, device_.getMaxWorkGroupSize());

    Buffer xBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*vecSize, const_cast<float*>(vecx));
    Buffer yBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*vecSize, const_cast<float*>(vecy));

    Buffer resultBuffer(context_, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(float));
    int32_t mutex = 0;
    Buffer mutexBuffer(context_, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(int32_t), &mutex);
    //Buffer mutexBuffer(context_, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(int32_t), &mutex);

    kernel->setArg(0, vecSize);
    kernel->setArg(1, xBuffer);
    kernel->setArg(2, yBuffer);
    kernel->setArg(3, resultBuffer);
    kernel->setArg(4, mutexBuffer);
    kernel->setArg(5, sizeof(float)*localWorkSize, 0);

    queue_->enqueue(kernel, workSize, localWorkSize);

    float result;
    queue_->enqueueReadBuffer(&resultBuffer, (void*)(&result), true);
    queue_->finish();

    return result;
}

float VoreenBlasCL::sNRM2(size_t vecSize, const float* vecx) const {

    if (!initialized_) {
        LERROR("Not initialized. Aborting.");
        return -1.f;
    }

    Kernel* kernel = prog_->getKernel("sNRM2");
    if (!kernel) {
        LERROR("No kernel 'sNRM2' found");
        return 0.f;
    }

    int workSize = 1024 << 2;
    int localWorkSize = std::min<int>(512, device_.getMaxWorkGroupSize());

    Buffer xBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*vecSize, const_cast<float*>(vecx));

    Buffer resultBuffer(context_, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float));
    int32_t mutex = 0;
    Buffer mutexBuffer(context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int32_t), &mutex);
    Buffer semaphorBuffer(context_, CL_MEM_ALLOC_HOST_PTR, sizeof(int32_t));

    kernel->setArg(0, vecSize);
    kernel->setArg(1, xBuffer);
    kernel->setArg(2, resultBuffer);
    kernel->setArg(3, mutexBuffer);
    kernel->setArg(4, semaphorBuffer);
    kernel->setArg(5, sizeof(float)*localWorkSize, 0);
    queue_->enqueue(kernel, workSize, localWorkSize);

    float result;
    queue_->enqueueReadBuffer(&resultBuffer, (void*)(&result), true);
    queue_->finish();

    return result;
}

void VoreenBlasCL::sSpMVEll(const EllpackMatrix<float>& mat, const float* vec, float* result) const {

    if (!initialized_) {
        LERROR("Not initialized. Aborting.");
        return;
    }

    Kernel* kernel = prog_->getKernel("sSpMV_Ell");
    if (!kernel) {
        LERROR("No kernel 'sSpMV_Ell' found");
        return;
    }

    int workSize = 1024 << 2;

    Buffer matBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, mat.getMatrixBufferSize(), mat.getMatrix());
    Buffer indicesBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, mat.getIndicesBufferSize(), mat.getIndices());
    Buffer vecBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*mat.getNumCols(), const_cast<float*>(vec));
    Buffer resultBuffer(context_, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*mat.getNumRows());

    kernel->setArg(0, mat.getNumRows());
    kernel->setArg(1, mat.getNumCols());
    kernel->setArg(2, mat.getNumColsPerRow());
    kernel->setArg(3, indicesBuffer);
    kernel->setArg(4, matBuffer);
    kernel->setArg(5, vecBuffer);
    kernel->setArg(6, resultBuffer);

    queue_->enqueue(kernel, workSize);

    queue_->enqueueReadBuffer(&resultBuffer, (void*)(result), true);
    queue_->finish();
}

void VoreenBlasCL::hSpMVEll(const EllpackMatrix<int16_t>& mat, const float* vec, float* result) const {

    if (!initialized_) {
        LERROR("Not initialized. Aborting.");
        return;
    }

    Kernel* kernel = prog_->getKernel("hSpMV_Ell");
    if (!kernel) {
        LERROR("No kernel 'hSpMV_Ell' found");
        return;
    }

    int workSize = 1024 << 2;

    Buffer matBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, mat.getMatrixBufferSize(), mat.getMatrix());
    Buffer indicesBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, mat.getIndicesBufferSize(), mat.getIndices());
    Buffer vecBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*mat.getNumCols(), const_cast<float*>(vec));
    Buffer resultBuffer(context_, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*mat.getNumRows());

    kernel->setArg(0, mat.getNumRows());
    kernel->setArg(1, mat.getNumCols());
    kernel->setArg(2, mat.getNumColsPerRow());
    kernel->setArg(3, indicesBuffer);
    kernel->setArg(4, matBuffer);
    kernel->setArg(5, vecBuffer);
    kernel->setArg(6, resultBuffer);

    queue_->enqueue(kernel, workSize);

    queue_->enqueueReadBuffer(&resultBuffer, (void*)(result), true);
    queue_->finish();
}

float VoreenBlasCL::sSpInnerProductEll(const EllpackMatrix<float>& mat, const float* vecx, const float* vecy) const {

    if (!initialized_) {
        LERROR("Not initialized. Aborting.");
        return -1.f;
    }

    Kernel* kernel = prog_->getKernel("sInnerProduct_Ell");
    if (!kernel) {
        LERROR("No kernel 'sInnerProduct_Ell' found");
        return -1.f;
    }

    int workSize = 1024;
    int localWorkSize = std::min<int>(256, device_.getMaxWorkGroupSize());

    Buffer matBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, mat.getMatrixBufferSize(), mat.getMatrix());
    Buffer indicesBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, mat.getIndicesBufferSize(), mat.getIndices());
    Buffer vecxBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*mat.getNumCols(), const_cast<float*>(vecx));
    Buffer vecyBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*mat.getNumCols(), const_cast<float*>(vecy));
    int32_t mutex = 0;
    Buffer mutexBuffer(context_, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(int32_t), &mutex);
    Buffer resultBuffer(context_, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(float));

    kernel->setArg(0, mat.getNumRows());
    kernel->setArg(1, mat.getNumCols());
    kernel->setArg(2, mat.getNumColsPerRow());
    kernel->setArg(3, indicesBuffer);
    kernel->setArg(4, matBuffer);
    kernel->setArg(5, vecxBuffer);
    kernel->setArg(6, vecyBuffer);
    kernel->setArg(7, resultBuffer);
    kernel->setArg(8, mutexBuffer);
    kernel->setArg(9, sizeof(float)*localWorkSize, 0);

    queue_->enqueue(kernel, workSize, localWorkSize);

    float result;
    queue_->enqueueReadBuffer(&resultBuffer, (void*)(&result), true);
    queue_->finish();

    return result;
}

int VoreenBlasCL::sSpConjGradEll(const EllpackMatrix<float>& mat, const float* vec, float* result,
                           float* initial, ConjGradPreconditioner precond, float threshold, int maxIterations, ProgressReporter& progress) const {

    if (!mat.isSymmetric()) {
        LERROR("Symmetric matrix expected.");
        return -1;
    }

    if (!initialized_) {
        LERROR("Not initialized. Aborting.");
        return -1;
    }

    Kernel* kernelSapxy = prog_->getKernel("sAXPY");
    if (!kernelSapxy) {
        LERROR("No kernel 'sAXPY' found");
        return -1;
    }

    Kernel* kernelSdot = prog_->getKernel("sDOT");
    if (!kernelSdot) {
        LERROR("No kernel 'sDOT' found");
        return -1;
    }

    Kernel* kernelSpmv = prog_->getKernel("sSpMV_Ell");
    if (!kernelSpmv) {
        LERROR("No kernel 'sSpMV_Ell' found");
        return -1;
    }

    Kernel* kernelSnrm2 = prog_->getKernel("sNRM2");
    if (!kernelSnrm2) {
        LERROR("No kernel 'sNRM2' found");
        return -1;
    }

    size_t vecSize = mat.getNumRows();
    int workSize = 1024 << 2;
    int localWorksize = std::min<int>(512, device_.getMaxWorkGroupSize());

    if (mat.getNumRows() < (size_t)device_.getMaxWorkGroupSize())
        workSize = localWorksize = device_.getMaxWorkGroupSize();

    bool initialAllocated = false;
    if (initial == 0) {
        initial = new float[vecSize];
        initialAllocated = true;
        for (size_t i=0; i<vecSize; ++i)
            initial[i] = 0.f;
    }

    EllpackMatrix<float>* preconditioner = 0;
    Buffer* precondBuf = 0;
    Buffer* precondIndicesBuf = 0;
    Buffer* zBuf = 0;

    Buffer matBuf(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, mat.getMatrixBufferSize(), mat.getMatrix());
    Buffer indicesBuf(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, mat.getIndicesBufferSize(), mat.getIndices());

    Buffer tmpBuf(context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*vecSize, const_cast<float*>(vec));
    Buffer xBuf(context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*vecSize, initial);
    Buffer rBuf(context_, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*vecSize);
    Buffer pBuf(context_, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*vecSize);

    int iteration  = 0;

    try {

        if (precond == Jacobi) {
            preconditioner = new EllpackMatrix<float>(mat.getNumRows(), mat.getNumRows(), 1);
            preconditioner->initializeBuffers();
            for (size_t i=0; i<mat.getNumRows(); i++)
                preconditioner->setValueByIndex(i, i, 0, 1.f / std::max(mat.getValue(i,i), 1e-6f));

            precondBuf = new Buffer(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                    preconditioner->getMatrixBufferSize(), preconditioner->getMatrix());
            precondIndicesBuf = new Buffer(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, preconditioner->getIndicesBufferSize(), preconditioner->getIndices());

            zBuf = new Buffer(context_, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*vecSize);
        }

        int32_t mutex;
        Buffer mutexBuf(context_, CL_MEM_ALLOC_HOST_PTR, sizeof(int32_t));
        Buffer scalarBuf(context_, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(float));
        Buffer nominatorBuf(context_, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float));
        Buffer denominatorBuf(context_, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float));
        Buffer alphaBuf(context_, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float));

        // r <= A*x_0
        kernelSpmv->setArg(0, vecSize);
        kernelSpmv->setArg(1, vecSize);
        kernelSpmv->setArg(2, mat.getNumColsPerRow());
        kernelSpmv->setArg(3, indicesBuf);
        kernelSpmv->setArg(4, matBuf);
        kernelSpmv->setArg(5, xBuf);
        kernelSpmv->setArg(6, rBuf);
        queue_->enqueue(kernelSpmv, workSize);

        // p <= -r + b
        kernelSapxy->setArg(0, vecSize);
        kernelSapxy->setArg(1, rBuf);
        kernelSapxy->setArg(2, tmpBuf);
        kernelSapxy->setArg(3, -1.f);
        kernelSapxy->setArg(4, pBuf);
        queue_->enqueue(kernelSapxy, workSize);

        // r <= -r + b
        kernelSapxy->setArg(0, vecSize);
        kernelSapxy->setArg(1, rBuf);
        kernelSapxy->setArg(2, tmpBuf);
        kernelSapxy->setArg(3, -1.f);
        kernelSapxy->setArg(4, rBuf);
        queue_->enqueue(kernelSapxy, workSize);

        // preconditioning
        if (precond == Jacobi) {
            kernelSpmv->setArg(0, vecSize);
            kernelSpmv->setArg(1, vecSize);
            kernelSpmv->setArg(2, preconditioner->getNumColsPerRow());
            kernelSpmv->setArg(3, precondIndicesBuf);
            kernelSpmv->setArg(4, precondBuf);
            kernelSpmv->setArg(5, rBuf);
            kernelSpmv->setArg(6, zBuf);
            queue_->enqueue(kernelSpmv, workSize);

            //memcpy(pBuf, zBuf, vecSize * sizeof(float));
            kernelSapxy->setArg(0, vecSize);
            kernelSapxy->setArg(1, zBuf);
            kernelSapxy->setArg(2, zBuf);
            kernelSapxy->setArg(3, 0.f);
            kernelSapxy->setArg(4, pBuf);
            queue_->enqueue(kernelSapxy, workSize);
        }

        while (iteration < maxIterations) {
            progress.setProgress(static_cast<float>(iteration)/maxIterations);

            iteration++;

            // r_k^T*r_k
            kernelSdot->setArg(0, vecSize);
            kernelSdot->setArg(1, rBuf);
            if (precond == Jacobi)
                kernelSdot->setArg(2, zBuf);
            else
                kernelSdot->setArg(2, rBuf);
            kernelSdot->setArg(3, nominatorBuf);
            kernelSdot->setArg(4, mutexBuf);
            kernelSdot->setArg(5, sizeof(float)*localWorksize, 0);
            mutex = 0;
            queue_->enqueueWriteBuffer(&mutexBuf, &mutex, true);    //< initialize mutex
            queue_->enqueue(kernelSdot, workSize, localWorksize);

            // tmp <= A * p_k
            kernelSpmv->setArg(0, vecSize);
            kernelSpmv->setArg(1, vecSize);
            kernelSpmv->setArg(2, mat.getNumColsPerRow());
            kernelSpmv->setArg(3, indicesBuf);
            kernelSpmv->setArg(4, matBuf);
            kernelSpmv->setArg(5, pBuf);
            kernelSpmv->setArg(6, tmpBuf);
            queue_->enqueue(kernelSpmv, workSize);

            // dot(p_k^T, tmp)
            kernelSdot->setArg(0, vecSize);
            kernelSdot->setArg(1, pBuf);
            kernelSdot->setArg(2, tmpBuf);
            kernelSdot->setArg(3, denominatorBuf);
            kernelSdot->setArg(4, mutexBuf);
            kernelSdot->setArg(5, sizeof(float)*localWorksize, 0);
            mutex = 0;
            queue_->enqueueWriteBuffer(&mutexBuf, &mutex, true);
            queue_->enqueue(kernelSdot, workSize, localWorksize);

            float denominator;
            float nominator;
            queue_->enqueueReadBuffer(&nominatorBuf, (void*)(&nominator), true);
            queue_->enqueueReadBuffer(&denominatorBuf, (void*)(&denominator), true);
            float alpha = nominator / denominator;

            // x <= alpha*p + x
            kernelSapxy->setArg(0, vecSize);
            kernelSapxy->setArg(1, pBuf);
            kernelSapxy->setArg(2, xBuf);
            kernelSapxy->setArg(3, alpha);
            kernelSapxy->setArg(4, xBuf);
            queue_->enqueue(kernelSapxy, workSize);

            // r <= -alpha*tmp + r
            kernelSapxy->setArg(0, vecSize);
            kernelSapxy->setArg(1, tmpBuf);
            kernelSapxy->setArg(2, rBuf);
            kernelSapxy->setArg(3, -alpha);
            kernelSapxy->setArg(4, rBuf);
            queue_->enqueue(kernelSapxy, workSize);

            // norm(r_k+1)
            float beta;
            if (precond == Jacobi) {
                kernelSpmv->setArg(0, vecSize);
                kernelSpmv->setArg(1, vecSize);
                kernelSpmv->setArg(2, preconditioner->getNumColsPerRow());
                kernelSpmv->setArg(3, precondIndicesBuf);
                kernelSpmv->setArg(4, precondBuf);
                kernelSpmv->setArg(5, rBuf);
                kernelSpmv->setArg(6, zBuf);
                queue_->enqueue(kernelSpmv, workSize);

                kernelSdot->setArg(0, vecSize);
                kernelSdot->setArg(1, rBuf);
                kernelSdot->setArg(2, zBuf);
                kernelSdot->setArg(3, scalarBuf);
                kernelSdot->setArg(4, mutexBuf);
                kernelSdot->setArg(5, sizeof(float)*localWorksize, 0);
                mutex = 0;
                queue_->enqueueWriteBuffer(&mutexBuf, &mutex, true);
                queue_->enqueue(kernelSdot, workSize, localWorksize);
                queue_->enqueueReadBuffer(&scalarBuf, (void*)(&beta), true);
            }
            else {
                kernelSdot->setArg(0, vecSize);
                kernelSdot->setArg(1, rBuf);
                kernelSdot->setArg(2, rBuf);
                kernelSdot->setArg(3, scalarBuf);
                kernelSdot->setArg(4, mutexBuf);
                kernelSdot->setArg(5, sizeof(float)*localWorksize, 0);
                mutex = 0;
                queue_->enqueueWriteBuffer(&mutexBuf, &mutex, true);
                queue_->enqueue(kernelSdot, workSize, localWorksize);
                queue_->enqueueReadBuffer(&scalarBuf, (void*)(&beta), true);
            }

            if (sqrt(beta) < threshold) {
                break;
            }

            beta /= nominator;

            // p <= beta*p + r
            kernelSapxy->setArg(0, vecSize);
            kernelSapxy->setArg(1, pBuf);
            if (precond == Jacobi)
                kernelSapxy->setArg(2, zBuf);
            else
                kernelSapxy->setArg(2, rBuf);
            kernelSapxy->setArg(3, beta);
            kernelSapxy->setArg(4, pBuf);
            queue_->enqueue(kernelSapxy, workSize);
        }

        queue_->enqueueReadBuffer(&xBuf, (void*)(result), true);

    } catch(boost::thread_interrupted& e) {
        queue_->finish();

        if (initialAllocated)
            delete[] initial;

        delete preconditioner;
        delete precondBuf;
        delete precondIndicesBuf;
        delete zBuf;

        throw e;
    }
    queue_->finish();
    progress.setProgress(1.0);

    if (initialAllocated)
        delete[] initial;

    delete preconditioner;
    delete precondBuf;
    delete precondIndicesBuf;
    delete zBuf;

    return iteration;
}

int VoreenBlasCL::hSpConjGradEll(const EllpackMatrix<int16_t>& mat, const float* vec, float* result,
                                   float* initial, float threshold, int maxIterations) const {

    if (!mat.isSymmetric()) {
        LERROR("Symmetric matrix expected.");
        return -1;
    }

    if (!initialized_) {
        LERROR("Not initialized. Aborting.");
        return -1;
    }

    Kernel* kernelSapxy = prog_->getKernel("sAXPY");
    if (!kernelSapxy) {
        LERROR("No kernel 'sAXPY' found");
        return -1;
    }

    Kernel* kernelSdot = prog_->getKernel("sDOT");
    if (!kernelSdot) {
        LERROR("No kernel 'sDOT' found");
        return -1;
    }

    Kernel* kernelSpmv = prog_->getKernel("hSpMV_Ell");
    if (!kernelSpmv) {
        LERROR("No kernel 'hSpMV_Ell' found");
        return -1;
    }

    size_t vecSize = mat.getNumRows();
    cl_uint workSize = 1024 << 2;
    int localWorksize = std::min<int>(512, device_.getMaxWorkGroupSize());

    if (mat.getNumRows() < device_.getMaxWorkGroupSize())
        workSize = localWorksize = device_.getMaxWorkGroupSize();

    bool initialAllocated = false;
    if (initial == 0) {
        initial = new float[vecSize];
        initialAllocated = true;
        for (size_t i=0; i<vecSize; ++i)
            initial[i] = 0.f;
    }

    Buffer matBuf(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, mat.getMatrixBufferSize(), mat.getMatrix());
    Buffer indicesBuf(context_, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, mat.getIndicesBufferSize(), mat.getIndices());

    Buffer tmpBuf(context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*vecSize, const_cast<float*>(vec));
    Buffer xBuf(context_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*vecSize, initial);
    Buffer rBuf(context_, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*vecSize);
    Buffer pBuf(context_, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*vecSize);

    uint32_t mutex;
    Buffer mutexBuf(context_, CL_MEM_ALLOC_HOST_PTR, sizeof(int32_t));
    Buffer scalarBuf(context_, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(float));
    Buffer nominatorBuf(context_, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float));
    Buffer denominatorBuf(context_, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float));
    Buffer alphaBuf(context_, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float));

    // r <= A*x_0
    kernelSpmv->setArg(0, vecSize);
    kernelSpmv->setArg(1, vecSize);
    kernelSpmv->setArg(2, mat.getNumColsPerRow());
    kernelSpmv->setArg(3, indicesBuf);
    kernelSpmv->setArg(4, matBuf);
    kernelSpmv->setArg(5, xBuf);
    kernelSpmv->setArg(6, rBuf);
    queue_->enqueue(kernelSpmv, workSize);

    // p <= -r + b
    kernelSapxy->setArg(0, vecSize);
    kernelSapxy->setArg(1, rBuf);
    kernelSapxy->setArg(2, tmpBuf);
    kernelSapxy->setArg(3, -1.f);
    kernelSapxy->setArg(4, pBuf);
    queue_->enqueue(kernelSapxy, workSize);

    // r <= -r + b
    kernelSapxy->setArg(0, vecSize);
    kernelSapxy->setArg(1, rBuf);
    kernelSapxy->setArg(2, tmpBuf);
    kernelSapxy->setArg(3, -1.f);
    kernelSapxy->setArg(4, rBuf);
    queue_->enqueue(kernelSapxy, workSize);

    int iteration  = 0;

    while (iteration < maxIterations) {

        iteration++;

        // r_k^T*r_k
        kernelSdot->setArg(0, vecSize);
        kernelSdot->setArg(1, rBuf);
        kernelSdot->setArg(2, rBuf);
        kernelSdot->setArg(3, nominatorBuf);
        kernelSdot->setArg(4, mutexBuf);
        kernelSdot->setArg(5, sizeof(float)*localWorksize, 0);
        mutex = 0;
        queue_->enqueueWriteBuffer(&mutexBuf, &mutex, true);    //< initialize mutex
        queue_->enqueue(kernelSdot, workSize, localWorksize);

        // tmp <= A * p_k
        kernelSpmv->setArg(0, vecSize);
        kernelSpmv->setArg(1, vecSize);
        kernelSpmv->setArg(2, mat.getNumColsPerRow());
        kernelSpmv->setArg(3, indicesBuf);
        kernelSpmv->setArg(4, matBuf);
        kernelSpmv->setArg(5, pBuf);
        kernelSpmv->setArg(6, tmpBuf);
        queue_->enqueue(kernelSpmv, workSize);

        // dot(p_k^T, tmp)
        kernelSdot->setArg(0, vecSize);
        kernelSdot->setArg(1, pBuf);
        kernelSdot->setArg(2, tmpBuf);
        kernelSdot->setArg(3, denominatorBuf);
        kernelSdot->setArg(4, mutexBuf);
        kernelSdot->setArg(5, sizeof(float)*localWorksize, 0);
        mutex = 0;
        queue_->enqueueWriteBuffer(&mutexBuf, &mutex, true);
        queue_->enqueue(kernelSdot, workSize, localWorksize);

        float denominator;
        float nominator;
        queue_->enqueueReadBuffer(&nominatorBuf, (void*)(&nominator), true);
        queue_->enqueueReadBuffer(&denominatorBuf, (void*)(&denominator), true);
        float alpha = nominator / denominator;

        // x <= alpha*p + x
        kernelSapxy->setArg(0, vecSize);
        kernelSapxy->setArg(1, pBuf);
        kernelSapxy->setArg(2, xBuf);
        kernelSapxy->setArg(3, alpha);
        kernelSapxy->setArg(4, xBuf);
        queue_->enqueue(kernelSapxy, workSize);

        // r <= -alpha*tmp + r
        kernelSapxy->setArg(0, vecSize);
        kernelSapxy->setArg(1, tmpBuf);
        kernelSapxy->setArg(2, rBuf);
        kernelSapxy->setArg(3, -alpha);
        kernelSapxy->setArg(4, rBuf);
        queue_->enqueue(kernelSapxy, workSize);

        // norm(r_k+1)
        kernelSdot->setArg(0, vecSize);
        kernelSdot->setArg(1, rBuf);
        kernelSdot->setArg(2, rBuf);
        kernelSdot->setArg(3, scalarBuf);
        kernelSdot->setArg(4, mutexBuf);
        kernelSdot->setArg(5, sizeof(float)*localWorksize, 0);
        mutex = 0;
        queue_->enqueueWriteBuffer(&mutexBuf, &mutex, true);
        queue_->enqueue(kernelSdot, workSize, localWorksize);
        float beta;
        queue_->enqueueReadBuffer(&scalarBuf, (void*)(&beta), true);

        if (sqrt(beta) < threshold)
            break;

        beta /= nominator;

        // p <= beta*p + r
        kernelSapxy->setArg(0, vecSize);
        kernelSapxy->setArg(1, pBuf);
        kernelSapxy->setArg(2, rBuf);
        kernelSapxy->setArg(3, beta);
        kernelSapxy->setArg(4, pBuf);
        queue_->enqueue(kernelSapxy, workSize);
    }

    queue_->enqueueReadBuffer(&xBuf, (void*)(result), true);
    queue_->finish();

    if (initialAllocated)
        delete[] initial;

    return iteration;
}

}   // namespace
