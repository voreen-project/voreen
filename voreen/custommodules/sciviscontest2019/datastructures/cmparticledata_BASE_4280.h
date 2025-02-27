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

#ifndef VRN_CMPARTICLEDATA_H
#define VRN_CMPARTICLEDATA_H

#include "tgt/vector.h"
#include "tgt/bounds.h"
#include <string>
#include <vector>

namespace voreen {

struct CMParticle {
    tgt::vec3 pos;	/* Positions (x, y, z) are in units of Mpc/h */
    tgt::vec3 vel;	/* Velocities (vx, vy, vz) are in units of km/s */
	float mass;		/* Mass is in internal code units */
	float uu;		/* Internal energy(uu) is in units of(km / s) ^ 2 */
	float hh;		/* SPH smoothing length(hh) is in units of Mpc / h */
	float mu;		/* Molecular weight (mu) is dimensionless */
	float rho;		/* Density (rho) is in units of h^2*Msolar/Mpc^3 where Msolar is the mass of the Sun */
    float phi;      /* Gravitational potential (phi) is in internal code units */
    int64_t ident;  /* unique identifier */
	uint16_t mask;	/*2nd bit : Denotes whether a particle is dark matter(0) or a baryon(1)
					  6th bit : Denotes if a baryon particle is also a star particle
					  7 bit : Denotes if a baryon particle is also a wind particle
					  8th bit : Denotes if a baryon particle is also a star - forming gas particle
					  9th bit : Denotes if a dark matter particle has been flagged as an active galactic nuclei(AGN) */
};

class CosmologyParticleFilter;
class CMParticleDataTimeSlice {
public:
    friend CosmologyParticleFilter;

    CMParticleDataTimeSlice(tgt::Bounds universeBounds);
    virtual const CMParticle* startUsingParticles()  = 0;
    virtual void              finishedUsingParticles() {};
    virtual int               getNumberOfParticles() = 0;
    virtual float             getTimeStep() = 0;
    virtual const int*        getRemappingTable() = 0;

    int*                      buildRemappingTable();


    virtual ~CMParticleDataTimeSlice(){}

    tgt::vec2 getMinMaxVelocityMagnitude();
    tgt::vec2 getMinMaxPhi();
    //After NormalisatoinTransformation => in Mpc/h
    virtual tgt::Bounds getBounds();
    //In original coordinates
    virtual tgt::Bounds getUniverseBounds();
    tgt::Bounds getVelocityBounds();
    virtual float geth0() = 0;
    //Scaling and shifting to comoving [0, L0]^3 in Mpc/h
    virtual tgt::mat4 getNormalisationTransformation();

    tgt::Bounds universeBounds_;
};

class CMParticleDataTimeSliceVector : public CMParticleDataTimeSlice{
public:
    CMParticleDataTimeSliceVector(tgt::Bounds universeBounds);
    virtual const CMParticle* startUsingParticles() override;
    virtual int               getNumberOfParticles() override;
    virtual float             getTimeStep() override;
    virtual const int*        getRemappingTable() override;
    virtual float             geth0() override;

    std::vector<CMParticle> particles_;

    float timeStep_;
    float h0_;

private:
    std::vector<int> remappingTable_;
};


class CMParticleData : public DataInvalidationObservable {
public:
    CMParticleData(const std::vector<CMParticleDataTimeSlice*>& timeSlices);
    ~CMParticleData();

    CMParticleData* cloneWithParticles() const;

    const std::vector<CMParticleDataTimeSlice*>& particleDataTimeSlices() const;
    CMParticleDataTimeSlice* sliceAtTimeStep(float a) const;

    const char* getEnabledState() const;
    char* getEnabledStateMutable();

    bool isFiltered() const;
private:
    std::vector<CMParticleDataTimeSlice *> timeSlices_;
    bool isFiltered_;
    std::vector<char> enabledState_;

    static const std::string loggerCat_; ///< category used in logging
};
} // namespace
#endif
