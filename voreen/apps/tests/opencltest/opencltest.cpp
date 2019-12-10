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

#include "GL/glew.h"
#include "tgt/logmanager.h"
#include "tgt/init.h"
#include "voreen/core/utils/clwrapper.h"

using namespace voreen;

int main() {
    tgt::init();
    const std::string loggerCat_ = "opencltest";

    OpenCL opencl;

    const std::vector<Platform>&  platforms = opencl.getPlatforms();
    if(platforms.size() == 0) {
        LERROR("Found no OpenCL platforms!");
        return -1;
    }

    const std::vector<Device*>& devices = platforms[0].getDevices();
    if(devices.size() == 0) {
        LERROR("Found no devices in platform!");
        return -1;
    }

    Context con(devices[0]);
    std::string source = "__kernel void sq(__global float* input, __global float* output, const unsigned int count) { int i = get_global_id(0); if(i < count) output[i] = input[i] * input[i]; }";

    Program prog(con);
    //prog.setSource(source);

    prog.loadSource("square.cl");

    //std::vector<std::string> filenames;
    //filenames.push_back("square.cl");
    //filenames.push_back("square2.cl");
    //prog.loadSource(filenames);

    prog.build(devices.back());
    Kernel* kernel = prog.getKernel("sq");
    if(kernel) {
        float in[20];
        float out[20];
        unsigned int count = 20;
        for(size_t i=0; i<count; ++i)
            in[i] = i;

        CommandQueue queue(&con, devices.back());

        //create buffers
        Buffer inBuffer(con, CL_MEM_READ_ONLY, sizeof(float) * count);
        Buffer outBuffer(con, CL_MEM_WRITE_ONLY, sizeof(float) * count);
        //upload input
        queue.enqueueWriteBuffer(&inBuffer, &in[0]);
        //set kernel args
        kernel->setArg(0, inBuffer);
        kernel->setArg(1, outBuffer);
        kernel->setArg(2, count);
        //run kernel
        queue.enqueue(kernel, count, 1);
        queue.finish();

        //read back results
        queue.enqueueReadBuffer(&outBuffer, &out[0]);
        queue.finish();
        //output results
        for(size_t i=0; i<count; ++i)
            LINFO("in[" << i << "]=" << in[i] << " out[" << i << "]=" << out[i]);
    }

    //prog.loadSource("square.cl");
    prog.setBuildOptions("-Dtest");
    prog.build(devices.back());
    kernel = prog.getKernel("sq");
    if(kernel) {
        float in[20];
        float out[20];
        unsigned int count = 20;
        for(size_t i=0; i<count; ++i)
            in[i] = i;

        CommandQueue queue(&con, devices.back());

        //create buffers
        Buffer inBuffer(con, CL_MEM_READ_ONLY, sizeof(in));
        Buffer outBuffer(con, CL_MEM_WRITE_ONLY, sizeof(out));
        //upload input
        queue.enqueueWriteBuffer(&inBuffer, in);
        //set kernel args
        kernel->setArg(0, inBuffer);
        kernel->setArg(1, outBuffer);
        kernel->setArg(2, count);
        //run kernel
        queue.enqueue(kernel, count, 1);
        queue.finish();

        //read back results
        queue.enqueueReadBuffer(&outBuffer, out);
        //output results
        for(size_t i=0; i<count; ++i)
            LINFO("in[" << i << "]=" << in[i] << " out[" << i << "]=" << out[i]);
    }
    return 0;
}
