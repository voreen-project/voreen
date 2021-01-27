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

#ifndef VRN_STREAMLINELISTDECORATOR_H
#define VRN_STREAMLINELISTDECORATOR_H

#include "streamlinelistbase.h"
#include "streamlinelistobserver.h"

namespace voreen {

/**
 * Decorator base class of SteamlineListDecorators.
 */
class VRN_CORE_API StreamlineListDecoratorIdentity : public StreamlineListBase, public StreamlineListObserver {
public:
    //------------------------------
    //  Constructor and Destructor
    //------------------------------

    /** Constructor */
    StreamlineListDecoratorIdentity(StreamlineListBase* base);
    /** Destructor */
    virtual ~StreamlineListDecoratorIdentity();
    /** No Copy-Constructor. Use clone instead. */
    virtual StreamlineListBase* clone() const;
private:
    /** No Copy */
    StreamlineListDecoratorIdentity(const StreamlineListDecoratorIdentity&);
    /** No Copy */
    StreamlineListDecoratorIdentity & operator=(const StreamlineListDecoratorIdentity&);
public:
    //------------------------------
    //  Observer Functions
    //------------------------------
    virtual void beforeStreamlineListDelete(const StreamlineListBase* source);

    //------------------------------
    //  Streamline Handling
    //------------------------------
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
    virtual const tgt::mat4   getVoxelToWorldMatrix() const;

    virtual float getMinMagnitude() const;
    virtual float getMaxMagnitude() const;

    virtual const tgt::vec2& getTemporalRange() const;

protected:
    virtual void setTransformMatrices(const tgt::mat4& listMatrix, const tgt::mat4& velocityMatrix);

public:
    //------------------------------
    //  Storage
    //------------------------------
    virtual std::string metaToCSVString() const;
    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    //------------------------------
    //  Members
    //------------------------------
protected:
    mutable StreamlineListBase* decorated_; ///< pointer to the base StreamlineListBase...
};




/**
 * Decorator to replace the transformation matrix from a SteamlineListBase.
 */
class VRN_CORE_API StreamlineListDecoratorReplaceTransformation : public StreamlineListDecoratorIdentity {
public:
    StreamlineListDecoratorReplaceTransformation(StreamlineListBase* base, const tgt::mat4& listMatrix, const tgt::mat4& velocityMatrix);
    virtual ~StreamlineListDecoratorReplaceTransformation(){}
     /** No Copy-Constructor. Use clone instead. */
    virtual StreamlineListBase* clone() const;
private:
    StreamlineListDecoratorReplaceTransformation(const StreamlineListDecoratorReplaceTransformation&);
    StreamlineListDecoratorReplaceTransformation & operator=(const StreamlineListDecoratorReplaceTransformation&);
public:
    virtual const tgt::mat4&   getListTransformMatrix() const;
    virtual const tgt::mat4&   getVelocityTransformMatrix() const;
    virtual const tgt::mat4   getVoxelToWorldMatrix() const;

    virtual void serialize(Serializer& s) const;
protected:
    virtual void setTransformMatrices(const tgt::mat4& listMatrix, const tgt::mat4& velocityMatrix);

    //------------------------------
    //  Members
    //------------------------------
    tgt::mat4 decoratorListTransformMatrix_;         ///< transforms the list/position
    tgt::mat4 decoratorVelocityTransformMatrix_;     ///< transforms each velocity
};

}   // namespace

#endif
