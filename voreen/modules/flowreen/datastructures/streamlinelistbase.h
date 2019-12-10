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

#ifndef VRN_STREAMLINELISTBASE_H
#define VRN_STREAMLINELISTBASE_H

#include "voreen/core/io/serialization/serializable.h"
#include "streamlinelistobserver.h"

#include "streamline.h"
#include "streamlinebundle.h"

#include "tgt/vector.h"
#include "tgt/matrix.h"
#include "tgt/bounds.h"

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API Observable<StreamlineListObserver>;
#endif

/**
 * A collection of Streamlines. It is used in the StreamlineListPort.
 * Beside the vector of Streamlines, meta-data of the original volume is stored. The streamlines themselves are
 * stored in the class StreamlineList
 * We use a base class, since a StreamlineListDecorator is used to override the transformation matrix.
 */
class VRN_CORE_API StreamlineListBase : public Serializable, public Observable<StreamlineListObserver> {
    friend class StreamlineRotation;
    friend class StreamlineListDecoratorIdentity;
    friend class StreamlineListDecoratorReplaceTransformation;
public:

    //------------------------------
    //  Constructor and Destructor
    //------------------------------

    /** Constructor */
    StreamlineListBase();
    /** Destructor */
    virtual ~StreamlineListBase();
    /** No Copy-Constructor. Use clone instead. */
    virtual StreamlineListBase* clone() const = 0;
private:
    /** No Copy */
    StreamlineListBase(const StreamlineListBase&);
    /** No Copy */
    StreamlineListBase & operator=(const StreamlineListBase&);

public:
    //------------------------------
    //  Observer Notification
    //------------------------------
    /** Notifies the registered StremlineListObservers about the pending deletion of the Volume. */
    void notifyDelete();

    //------------------------------
    //  Streamline Handling
    //------------------------------
    /** Adds a Streamline and copies it. */
    virtual void addStreamline(const Streamline& line) = 0;
    /** Adds and copies all Streamlines. No meta data is copied except min/max magnitude. */
    virtual void addStreamlineList(const StreamlineListBase& list) = 0;
    /** Removes a Streamline and returns all remaining ones. */
    virtual const std::vector<Streamline>& removeStreamline(size_t pos) = 0;
    /** Returns all Streamlines. */
    virtual const std::vector<Streamline>& getStreamlines() const = 0;

    //----------------------------
    //  Streamline Bundle Handling
    //----------------------------
    /** Adds a Streamline Bundle and copies it. */
    virtual void addStreamlineBundle(const StreamlineBundle& bundle) = 0;
    /** Removes a Streamline Bundle and returns all remaining ones. */
    virtual const std::vector<StreamlineBundle>& removeStreamlineBundle(size_t pos) = 0;
    /** Returns all Streamline Bundles. */
    virtual const std::vector<StreamlineBundle>& getStreamlineBundles() const = 0;
    /** Classifies a given streamline in terms of being noise in relation to bundles. */
    virtual void setStreamlineNoiseFlag(size_t pos) = 0;
    /** Returns all Streamlines being classified as noise. */
    virtual const std::vector<size_t>& getStreamlineNoise() const = 0;

    //------------------------------
    //  Meta
    //------------------------------
    virtual const tgt::svec3&  getOriginalDimensions() const = 0;
    virtual const tgt::vec3&   getOriginalSpacing() const  = 0;
    virtual const tgt::Bounds& getOriginalWorldBounds() const = 0;
    virtual const tgt::Bounds& getOriginalVoxelBounds() const = 0;
    virtual const tgt::mat4&   getOriginalVoxelToWorldMatrix() const = 0;
    virtual const tgt::mat4&   getOriginalWorldToVoxelMatrix() const = 0;

    virtual const tgt::mat4&   getListTransformMatrix() const = 0;
    virtual const tgt::mat4&   getVelocityTransformMatrix() const = 0;
    virtual const tgt::mat4    getVoxelToWorldMatrix() const = 0;

    virtual float getMinMagnitude() const = 0;
    virtual float getMaxMagnitude() const = 0;

protected:
    /**
     * Only used by the StreamlineRotation processor.
     * @param listMatrix Matrix for rotating the entire List. This includes translations!
     * @param velocityMatrix Matrix for the rotation of each element.
     */
    virtual void setTransformMatrices(const tgt::mat4& listMatrix, const tgt::mat4& velocityMatrix) = 0;

public:
    //------------------------------
    //  Storage
    //------------------------------
    virtual std::string metaToCSVString() const = 0;
    virtual void serialize(Serializer& s) const override = 0;
    virtual void deserialize(Deserializer& s) override = 0;
};

}   // namespace

#endif
