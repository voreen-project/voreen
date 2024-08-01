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

#ifndef VRN_CMMERGERTREE_H
#define VRN_CMMERGERTREE_H

#include <iostream>
#include <fstream>
#include "tgt/vector.h"
#include <vector>

namespace voreen{

class CMHalo;

/**
 * Represents a collection of halos in the same time step of a dark matter universe simulation.
 * @note This implementation relies heavily on the relative memory positions of objects of this
 *          and other classes in this file. This is ensured when the CMMergerTree is created by
 *          CMHaloDataSource.
 */
class CMTimeStepHaloList {
friend class CMHaloDataSource;
public:

    /// Actual time of this time step
    float a;

    /**
     * Begin of the memory block of halos that belong to this time step.
     */
    const CMHalo* halosBegin() const;

    /**
     * End of the memory block of halos that belong to this time step.
     * Points to the Halo just behind the last one that belongs to
     * this time step.
     */
    const CMHalo* halosEnd() const;

    /**
     * Returns the number of halos in this time step.
     */
    size_t size() const;

private:
    /**
     * Disable copy constructors because they break access to halos via offset:
     */
    CMTimeStepHaloList(const CMTimeStepHaloList&);

    /// Byte offset between this object and the first halo of this time step.
    size_t beginOffset;
    /// Byte offset between this object and the first halo NOT part of this time step:
    size_t endOffset;
};

/**
 * Represents a collection of halos of a dark matter universe simulation.
 * @note This implementation relies heavily on the relative memory positions of objects of this
 *          and other classes in this file. This is ensured when the CMMergerTree is created by
 *          CMHaloDataSource.
 */
class CMMergerTree
{
friend class CMHaloDataSource;
public:
    /// Symbolic ID that indicates that there is no such halo
    static const int NO_HALO_ID;

    /**
     * Constructor.
     */
    CMMergerTree();
    /**
     * Destructor.
     */
    ~CMMergerTree();

    /**
     * Returns the halos that exist in the time step closest to a.
     */
    const CMTimeStepHaloList* haloDataAt(float a) const;

    /**
     * Returns the halo identifiable by the given id. If there is
     * no such halo nullptr will be returned.
     */
    const CMHalo* haloByID(int id) const;

    /**
     * Begin of the memory block of time steps of this simulation.
     */
    const CMTimeStepHaloList* stepListsBegin() const;

    /**
     * End of the memory block of time steps that belong to this simulation.
     * Points to the time step just behind the last one that belongs
     * to this simulation.
     */
    const CMTimeStepHaloList* stepListsEnd() const;

    /**
     * Begin of the memory block of halos of this simulation.
     */
    const CMHalo* halosBegin() const;

    /**
     * End of the memory block of halos that belong to this simulation.
     * Points to the halo just behind the last one that belongs
     * to this simulation.
     */
    const CMHalo* halosEnd() const;

    /**
     * Indicates whether or not the halos belonging to this merger tree
     * contain additional data calculated by galacticus or not.
     * If false, the halos of course still have the members, but contain
     * arbitrary values.
     * @return Do halos contain additional attributes calculated by galacticus?
     */
    bool containsGalacticusData() const;

private:
    /**
     * Disable copy constructors because they break access to halos via offset.
     */
    CMMergerTree(const CMMergerTree&);

    /// Byte offset between this object and the first time step of the simulation.
    size_t hlistBeginOffset;
    /// Byte offset between this object and the first time step NOT part of this simulation.
    size_t hlistEndOffset;
    /// Byte offset between this object and the first halo of the simulation.
    size_t halosBeginOffset;
    /// Byte offset between this object and the first halo NOT part of this simulation.
    size_t halosEndOffset;
    /// Wether or not galacticus calculated attributes of halos contain meaningful values or not.
    bool   containsGalacticusData_;
};

/**
 * Represents a single halos of a dark matter universe simulation.
 * @note This implementation relies heavily on the relative memory positions of objects of this
 *          and other classes in this file. This is ensured when the CMMergerTree is created by
 *          CMHaloDataSource.
 */
class CMHalo
{
/// Let CMHaloDataSource create halo objects
friend class CMHaloDataSource;

public:
    
    /// The "time step". Actually probably the expansion factor that is related to, but not
    /// the same as the actual time. I'm not a physicist though.
    float scale;
    /// The unique identifier of this halo
    int ID;
    /// The unique identifier assigned by the simulation.
    int origID;
    /// The "time step" of the descendant of this halo
    float descendantScale;
    /// The identifier of the descendant of this halo
    int descendantID;
    /// The identifier of the host of this halo or NO_HALO_ID if this halo is not a satellite
    int hostID;
    /// The identifier of the (host's)* host of this halo or NO_HALO_ID if this halo is not a satellite
    int rootHostID;
    /// Original ID of the last descendant of this halo (in the last time step)
    int origTreeRootID;
    /// Mass of this halo.
    float mass;
    /// Radius of this halo.
    float radius;
    /// Scale radius of this halo.
    float scaleRadius;
    /// Position of this halo
    tgt::vec3 pos;
    /// Velocity of this halo
    tgt::vec3 velocity;
    /// Angular momentum of this halo
    tgt::vec3 angularMomenta;
    /// Spin parameter of this halo
    float spinParameter;

    /// ID of the first predecessor of this halo
    int parentID;
    /// ID of the next predecessor of this halo's descendant
    int spouseID;

    /// ID of the first satellite of this halo
    int satelliteID;
    /// ID of the next satellite of this halo's host
    int siblingSatelliteID;

    /// All calculated by galacticus
    double blackHoleMass;
    double blackHoleSpin;
    double spheroidRadius;
    double spheroidMassGas;
    double spheroidVelocity;
    double diskRadius;
    double diskMassGas;
    double diskVelocity;

    /// Pointer to this halo's descendant
    const CMHalo* descendant() const;
    /// Pointer to this halo's host or nullptr if there is no such halo
    const CMHalo* host() const;
    /// Pointer to this halo's (host's) host or nullptr if there is no such halo
    const CMHalo* rootHost() const;
    /// Pointer to this halo's first predecessor
    const CMHalo* parent() const;
    /// Pointer to this halo's descendant's next predecessor
    const CMHalo* spouse() const;
    /// Pointer to this halo's first satellite
    const CMHalo* satellite() const;
    /// Pointer to this halo's root's next satellite
    const CMHalo* siblingSatellite() const;

private:
    //CMHalo(const CMHalo&);
    CMHalo() {}
    /**
     * Get the address of the halo with the given id. No range checks will be applied.
     * @param id ID of the returned halo if valid!
     * @returns the halo identified by id or a pointer to somewhere else if id is not valid.
     */
    const CMHalo* unsafehaloByID(int id) const;
};
}
#endif //VRN_CMMERGERTREE_H
