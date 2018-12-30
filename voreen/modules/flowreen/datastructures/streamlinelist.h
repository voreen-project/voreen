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
//    friend class StreamlineRotation;
public:

    //------------------------------
    //  Constructor and Destructor
    //------------------------------

    /** Constructor, gets the meta informations from the passed volume */
    StreamlineList(const VolumeBase* vol = 0);
    /** Destructor */
    virtual ~StreamlineList();
    /** No Copy-Constructor. Use clone instead. */
    virtual StreamlineListBase* clone() const override;
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
    virtual void addStreamline(const Streamline& line) override;
    /** Adds and copies all Streamlines. No meta data is copied except min/max magnitude. */
    virtual void addStreamlineList(const StreamlineListBase& list) override;
    /** Removes a Streamline and returns all remaining ones. */
    virtual const std::vector<Streamline>& removeStreamline(size_t pos) override;
    /** Returns all Streamlines. */
    virtual const std::vector<Streamline>& getStreamlines() const override;

    //----------------------------
    //  Streamline Bundle Handling
    //----------------------------
    /** Adds a Streamline Bundle and copies it. */
    void addStreamlineBundle(const StreamlineBundle& bundle) override;
    /** Removes a Streamline Bundle and returns all remaining ones. */
    const std::vector<StreamlineBundle>& removeStreamlineBundle(size_t pos) override;
    /** Returns all Streamline Bundles. */
    const std::vector<StreamlineBundle>& getStreamlineBundles() const override;
    /** Classifies a given streamline in terms of being noise in relation to bundles. */
    void setStreamlineNoiseFlag(size_t pos) override;
    /** Returns all Streamlines being classified as noise. */
    const std::vector<size_t>& getStreamlineNoise() const override;

    //------------------------------
    //  Meta
    //------------------------------
    virtual const tgt::svec3&  getOriginalDimensions() const override;
    virtual const tgt::vec3&   getOriginalSpacing() const override;
    virtual const tgt::Bounds& getOriginalWorldBounds() const override;
    virtual const tgt::Bounds& getOriginalVoxelBounds() const override;
    virtual const tgt::mat4&   getOriginalVoxelToWorldMatrix() const override;
    virtual const tgt::mat4&   getOriginalWorldToVoxelMatrix() const override;

    virtual const tgt::mat4&   getListTransformMatrix() const override;
    virtual const tgt::mat4&   getVelocityTransformMatrix() const override;
    virtual const tgt::mat4    getVoxelToWorldMatrix() const override;

    virtual float getMinMagnitude() const override;
    virtual float getMaxMagnitude() const override;

protected:
    /**
     * Only used by the StreamlineRotation processor.
     * @param listMatrix Matrix for rotating the entire List. This includes translations!
     * @param velocityMatrix Matrix for the rotation of each element.
     */
    virtual void setTransformMatrices(const tgt::mat4& listMatrix, const tgt::mat4& velocityMatrix) override;

public:
    //----------------
    //  Storage
    //----------------
    virtual std::string metaToCSVString() const override;
    virtual void serialize(Serializer& s) const override;
    virtual void deserialize(Deserializer& s) override;

    //----------------
    //  Members
    //----------------
protected:

    std::vector<Streamline> streamlines_;               ///< list of streamlines
    std::vector<StreamlineBundle> streamlineBundles_;   ///< list of streamlinebundles
    std::vector<size_t> streamlineNoise_;               ///< list of indices classifying streamlines as noise

    //meta
    tgt::svec3 dimensions_;               ///< original dimensions of the input volume
    tgt::vec3 spacing_;                   ///< original spacing of the input volume
    tgt::Bounds worldBounds_;             ///< original bounds of the input volume (in world space)
    tgt::Bounds voxelBounds_;             ///< original bounds in voxel space
    tgt::mat4 voxelToWorldMatrix_;        ///< voxel to world matix
    tgt::mat4 worldToVoxelMatrix_;        ///< world to voxel matix
    float minMagnitude_;                  ///< global min magnitude value for color map configuration
    float maxMagnitude_;                  ///< global max magnitude value for color map configuration

    //Only used by the StreamlineRotation processor
    tgt::mat4 listTransformMatrix_;         ///< transforms the list/position
    tgt::mat4 velocityTransformMatrix_;     ///< transforms each velocity
};

}   // namespace

#endif
