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

#ifndef VRN_PREINTEGRATIONTABLE_H
#define VRN_PREINTEGRATIONTABLE_H

#include "tgt/vector.h"
#include "tgt/texture.h"
#include "voreen/core/processors/imageprocessor.h"
#include "voreen/core/utils/glsl.h"

#include "tgt/textureunit.h"
namespace voreen {

    class TransFunc1D;
/**
 * This realizes a pre-integration table (pre-integrated transfer function) as described by Engel et al., 2001.
 */
class VRN_CORE_API PreIntegrationTable
{
public:

    /**
     * Constructor, automatically calls computeTable at the end.
     *
     * @param transFunc the transfer function for with the pre-integration table is computed
     * @param resolution the resolution of the table in both directions (if <= 1 the resolution is set to 256)
     * @param d the segment length (= sampling step size) for which the pre-integration table is computed (if <= 0 the segment length is set to 1.0)
     * @param useIntegral Use integral functions to compute pre-integration table, which is faster but not quite as accurate.
     * @param computeOnGPU Compute the pre-integration table texture on the GPU.
     * @param program the shader program that is used to compute the pre-integration table on the GPU.
     */
    PreIntegrationTable(TransFunc1D* transFunc, size_t resolution = 256, float d = 1.f, bool useIntegral = true, bool computeOnGPU = false, tgt::Shader* program = 0);

    virtual ~PreIntegrationTable();

    /**
     * Compute the pre-integrated table for the given transfer function.
     */
    void computeTable() const;

    /**
     * Compute the pre-integrated table on the GPU.
     */
    void computeTableGPU() const;

    /**
     * Return the value in the pre-integration table.
     * This will return a vec4(0, 0, 0, 0) if the access is out of range [0, 1] or if no table is found.
     *
     * @param fs normalized intensity value at the start of the segment
     * @param fe normalized intensity value at the end of the segment
     */
    tgt::vec4 classify(float fs, float fe) const;

    /**
     * Returns the pre-integration table. If the tyble has not been computed yet, it is created (lazy computation).
     * For pre-integration tables computed on the GPU, this method returns 0. In that case, download the texture computed on the GPU
     */
    const tgt::vec4 * const getTable() const;

    /**
     * Returns the texture of the pre-integration table.
     * If no texture has been created yet, createTex() is called.
     */
    const tgt::Texture* getTexture() const;

    float getSamplingStepSize() const;

    bool usesIntegral() const;

    size_t getDimension() const;

private:

    // no default constructor
    PreIntegrationTable();

    /// Generates the texture for the pre-integration table.
    virtual void createTexFromTable() const;

    TransFunc1D* transFunc_;     ///< the 1D transfer function that is used to compute the pre-integration table

    size_t resolution_;         ///< resolution of the pre-integrated table
    float samplingStepSize_;    ///< length of the segments
    bool useIntegral_;          ///< true for approximative (but faster) computation using integral functions
    bool computeOnGPU_;         ///< if true this pre-integration table is computed on the GPU

    mutable tgt::vec4* table_; ///< the actual pre-integration table in row-major order

    mutable tgt::Texture* tex_; ///< texture for the pre-integration table, is generated internally

    mutable tgt::Shader* program_; ///< shader program to compute the pre-integration table on the gpu
    mutable RenderTarget renderTarget_; ///< internal render target for computing the pre-integration table on the gpu
};

} //namespace

#endif // VRN_PREINTEGRATIONTABLE_H
