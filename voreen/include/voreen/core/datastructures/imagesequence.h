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

#ifndef VRN_IMAGESEQUENCE_H
#define VRN_IMAGESEQUENCE_H

#include "voreen/core/voreencoreapi.h"

#include <vector>
#include <string>

namespace tgt {
    class Texture;
}

namespace voreen {

class VRN_CORE_API ImageSequence {

public:

    ImageSequence();

    virtual ~ImageSequence();

    /**
     * Adds the passed texture to the sequence, allowing duplicates.
     *
     * @note The passed texture must be a 2D texture.
     *
     * The ImageSequence does not take ownership of the
     * added texture and does therefore not delete it
     * on its own destruction.
     */
    virtual void add(tgt::Texture* texture);

    /**
     * Removes the passed texture from the sequence
     * without deleting it.
     */
    virtual void remove(const tgt::Texture* texture);

    /**
     * Removes the texture at position \p i from the sequence
     * without deleting it. \p i must be valid, i.e., i.e. i < size()
     */
    virtual void remove(size_t i);

    /**
     * Returns whether the passed texture is contained
     * by the ImageSequence.
     */
    virtual bool contains(const tgt::Texture* texture) const;

    /**
     * Returns the texture at a specified index position.
     *
     * @param i the index of the texture to return. Must be valid, i.e. i < size().
     */
    virtual tgt::Texture* at(size_t i) const;

    /**
     * Returns the first texture of the sequence, or null
     * if it is empty.
     */
    virtual tgt::Texture* front() const;

    /**
    * Returns the last element of the ImageSequence, or null
    * if it is empty.
    */
    virtual tgt::Texture* back() const;

    /**
     * Clears the ImageSequence without deleting the textures.
     */
    virtual void clear();

    /**
     * Returns the number of textures contained by the sequence,
     * including duplicates.
     */
    virtual size_t size() const;

    /**
     * Returns whether the ImageSequence is empty.
     */
    virtual bool empty() const;

protected:

    /// Vector storing the textures contained by the sequence.
    std::vector<tgt::Texture*> textures_;

    /// category for logging.
    static const std::string loggerCat_;

};

}   // namespace

#endif
