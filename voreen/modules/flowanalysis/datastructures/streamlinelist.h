/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_STREAMLINELIST_H
#define VRN_STREAMLINELIST_H

#include "streamlinelistbase.h"

#include "streamline.h"
#include "streamlinebundle.h"

#include <vector>
#include <string>

namespace voreen {

class VolumeBase;

/**
 * Implementation of StreamlineListBase, storing the streamlines.
 */
class VRN_CORE_API StreamlineList : public StreamlineListBase {
public:

    /** Constructor, gets the meta informations from the passed volume */
    StreamlineList(const VolumeBase* volume = nullptr);
    /** Destructor */
    virtual ~StreamlineList();
    /** No Copy-Constructor. Use clone instead. */
    virtual StreamlineListBase* clone() const;
private:
    /** No Copy */
    StreamlineList(const StreamlineList&);
    /** No Copy */
    StreamlineList & operator=(const StreamlineList&);

public:

    //------------------------
    //  Streamline Handling
    //------------------------
    /** Adds a Streamline and copies it. */
    virtual void addStreamline(const Streamline& line);
    /** Adds and copies all Streamlines. No meta data is copied except min/max magnitude. */
    virtual void addStreamlineList(const StreamlineListBase& list);
    /** Removes a Streamline and returns all remaining ones. */
    virtual void removeStreamline(size_t pos);
    /** Removes all streamlines */
    virtual void clearStreamlines() ;
    /** Returns all Streamlines. */
    virtual const std::vector<Streamline>& getStreamlines() const;

    //------------------------------
    //  Meta
    //------------------------------
    virtual const tgt::svec3&  getOriginalDimensions() const;
    virtual const tgt::vec3&   getOriginalSpacing() const;
    virtual const tgt::Bounds& getOriginalWorldBounds() const;
    virtual const tgt::Bounds& getOriginalVoxelBounds() const;
    virtual const tgt::mat4&   getOriginalVoxelToWorldMatrix() const;
    virtual const tgt::mat4&   getOriginalWorldToVoxelMatrix() const;

    virtual const tgt::mat4&   getListTransformMatrix() const;
    virtual const tgt::mat4&   getVelocityTransformMatrix() const;
    virtual const tgt::mat4    getVoxelToWorldMatrix() const;

    virtual float getMinMagnitude() const;
    virtual float getMaxMagnitude() const;

    virtual const tgt::vec2& getTemporalRange() const;

protected:
    /**
     * Only used by the StreamlineRotation processor.
     * @param listMatrix Matrix for rotating the entire List. This includes translations!
     * @param velocityMatrix Matrix for the rotation of each element.
     */
    virtual void setTransformMatrices(const tgt::mat4& listMatrix, const tgt::mat4& velocityMatrix);

public:
    //----------------
    //  Storage
    //----------------
    virtual std::string metaToCSVString() const;
    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    //----------------
    //  Members
    //----------------
protected:

    std::vector<Streamline> streamlines_; ///< list of streamlines

    //meta
    tgt::svec3 dimensions_;             ///< original dimensions of the input volume
    tgt::vec3 spacing_;                 ///< original spacing of the input volume
    tgt::Bounds worldBounds_;           ///< original bounds of the input volume (in world space)
    tgt::Bounds voxelBounds_;           ///< original bounds in voxel space
    tgt::mat4 voxelToWorldMatrix_;      ///< voxel to world matix
    tgt::mat4 worldToVoxelMatrix_;      ///< world to voxel matix

    // Statistics for rendering/color map adjustments
    tgt::vec2 magnitudeRange_;          ///< global magnitude range
    tgt::vec2 temporalRange_;           ///< global temporal range

    //Only used by the StreamlineRotation processor
    tgt::mat4 listTransformMatrix_;         ///< transforms the list/position
    tgt::mat4 velocityTransformMatrix_;     ///< transforms each velocity
};

}   // namespace

#endif
