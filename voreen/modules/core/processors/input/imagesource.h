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

#ifndef VRN_IMAGESOURCE_H
#define VRN_IMAGESOURCE_H

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/vectorproperty.h"

using tgt::Texture;

namespace voreen {

/**
 * Loads an image file from disk and supplies it as rendering.
 */
class VRN_CORE_API ImageSource : public RenderProcessor {

public:
    ImageSource();
    ~ImageSource();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "ImageSource";     }
    virtual std::string getCategory() const  { return "Input";           }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE; }

    /**
     * Loads the image specified by filename.
     *
     * @note The processor must have been initialized
     *       before calling this function.
     *
     * @param filename the image to load
     */
    void loadImage(const std::string& filename);

    /**
     * Assigns an image to this processor.
     *
     * @note The processor must have been initialized
     *       before calling this function.
     *
     * @param texture the 2D texture encapsulating the image.
     */
    void setImage(tgt::Texture* texture);

    /**
     * Returns the currently assigned image texture.
     */
    const tgt::Texture* getImage() const;

protected:
    virtual void setDescriptions() {
        setDescription("Loads an image file from disk and provides it as rendering. The outport dimensions are determined by the image dimensions.");
    }

    virtual void beforeProcess();
    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    /**
     * Clears the loaded image. Processor needs to be initialized when calling this function.
     */
    void clearImage();

    /**
     * Forces a reload of the loaded image.
     */
    void forceReload();

    RenderPort outport_;            ///< The render port the image is written to.
    FileDialogProperty imageFile_;  ///< Path of the loaded image file.
    ButtonProperty reloadImage_;    ///< Forces reload
    ButtonProperty clearImage_;     ///< Executes clearImage().
    IntVec2Property imageSize_;     ///< Read-only property providing the dimensions of the loaded image.

    tgt::Texture* texture_;         ///< The 2D texture encapsulating the image file.
    tgt::Shader* shader_;           ///< Shader program used for rendering the loaded image to the outport.

    /// Determines whether the processor owns the texture, which
    /// is only the case, if it has loaded it.
    bool textureOwner_;

    bool forceReload_;

    static const std::string loggerCat_;
};

} // namespace

#endif
