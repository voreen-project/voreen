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

#ifndef VRN_IMAGESEQUENCESOURCE_H
#define VRN_IMAGESEQUENCESOURCE_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/datastructures/imagesequence.h"

namespace voreen {

/**
 * Loads all image files from a directory and puts them out as ImageSequence
 * containing one OpenGL texture per image.
 */
class VRN_CORE_API ImageSequenceSource : public RenderProcessor {

public:
    ImageSequenceSource();
    ~ImageSequenceSource();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "ImageSequenceSource"; }
    virtual std::string getCategory() const  { return "Input";               }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE;     }

    /**
     * Issues a reload of the sequence during next beforeProcess(), if the directory path has changed.
     */
    virtual void invalidate(int inv = INVALID_RESULT);

    /**
     * Loads all images contained by the passed directory.
     *
     * @note The processor must have been initialized
     *       before calling this function.
     *
     * @param dir the directory containing the images to load
     */
    virtual void loadImageSequence(const std::string& dir);

    /**
     * Assigns the passed ImageSequence to this processor.
     *
     * @note The processor must have been initialized
     *       before calling this function.
     */
    virtual void setImageSequence(ImageSequence* sequence);

    /**
     * Returns the currently assigned/hold image sequence.
     */
    const ImageSequence* getImageSequence() const;

    /**
     * Issues a reload of the image sequence during next beforeProcess() call.
     */
    virtual void reloadImageSequence();

protected:
    virtual void setDescriptions() {
        setDescription("Loads all images from the specified directory and puts them out as image sequence containing one OpenGL texture per image.\
<p>See ImageSelector.</p>");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    /**
     * Clears the assigned image sequence.
     *
     * @note The processor must have been initialized
     *       before calling this function.
     */
    virtual void clearSequence();

    /**
     * Enforces a reload of the current image sequence during next beforeProcess() call.
     */
    virtual void forceReload();

    /**
     * Sets the imageDirectory_ property to an empty string,
     * thereby causing the sequence to be cleared during next beforeProcess().
     */
    virtual void unsetDirectoryName();

    ImageSequencePort outport_;         ///< The port the generated image sequence is written to.

    FileDialogProperty imageDirectory_; ///< Directory the images are loaded from.
    BoolProperty textureFiltering_;     ///< Enable linear filtering of loaded textures.
    BoolProperty showProgressBar_;      ///< Determines whether a progress dialog is shown during image loading.
    ButtonProperty reloadSequence_;     ///< Button for reloading the current sequence.
    ButtonProperty clearSequence_;      ///< Button for clearing the current sequence.
    IntProperty numImages_;             ///< Read-only property providing the size of the sequence.

    ImageSequence* imageSequence_;      ///< Image sequence assigned to this processor.
    std::string currentDir_;            ///< Current image directory.

    /// Determines whether the processor owns the sequence's textures,
    /// which is only the case, if it has loaded them.
    bool sequenceOwner_;

    bool forceReload_;
    static const std::string loggerCat_;
};

} // namespace

#endif
