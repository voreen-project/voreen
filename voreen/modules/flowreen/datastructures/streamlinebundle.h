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

#ifndef VRN_STREAMLINEBUNDLE_H
#define VRN_STREAMLINEBUNDLE_H

#include "streamline.h"

#include "voreen/core/utils/exception.h"
#include "voreen/core/voreencoreapi.h"
#include "voreen/core/io/serialization/serializable.h"

#include "tgt/vector.h"
#include "tgt/matrix.h"

#include <vector>

namespace voreen {

    /**
     * This class represents a bundle of streamlines, therefore storing the index of every streamline
     * being assigned. Those indices refer to the StreamlineList, the bundles are being created on.
     */
    class VRN_CORE_API StreamlineBundle : public Serializable {
    public:

        /**
         * Within the current implementation, a node contains no more information than a StreamlineElement.
         * Thus, we use the existing code and add information later on, if needed.
         */
        typedef Streamline::StreamlineElement Node;

    public:

        StreamlineBundle();
        StreamlineBundle(size_t index, const Streamline& prototype);

        /** Assignes a streamline to this bundle. */
        void addStreamline(size_t index, const Streamline& streamline);

        /** Returns the indices of all assigned streamlines. */
        const std::vector<size_t>& getStreamlines() const;

        /** Returns the Centroid Streamline of this bundle. */
        Streamline getCentroid() const;

        /** Returns a radius representing the extrution of the centroid. */
        float getRadius() const;

        //----------------
        //  Storage
        //----------------
        /** Used to save as CSV file. Transforms the voxel position to world space. */
        std::string toCSVString(const tgt::mat4& transfomationMatrix = tgt::mat4::identity) const;
        /** @override */
        virtual void serialize(Serializer& s) const;
        /** @override */
        virtual void deserialize(Deserializer& s);

    private:

        std::vector<size_t> streamlines_; ///< list containing all assigned streamlines
        std::vector<Node> nodes_;         ///< list containing all nodes
        float radius_;                    ///< a radius representing the extrution of the centroid

    };

}

#endif // VRN_STREAMLINEBUNDLE_H
