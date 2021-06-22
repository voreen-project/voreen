/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Mathias J. Krause, Marie-Luise Maier
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#include "dynamics/latticeDescriptors.h"

#include "functors/lattice/latticeInterpPhysVelocity3D.h"
#include "functors/lattice/latticeInterpPhysVelocity3D.hh"
#include "particle3D.h"
#include "particle3D.hh"
#include "superParticleSystem3D.h"
#include "superParticleSystem3D.hh"

namespace olb {

template class SuperParticleSystem3D<double, Particle3D>;

template class SuperParticleSystem3D<double, MagneticParticle3D>;

template<>
template<>
void SuperParticleSystem3D<double, Particle3D>::
setVelToFluidVel<descriptors::D3Q19<>>(
                                      SuperLatticeInterpPhysVelocity3D<double, descriptors::D3Q19<>>& fVel)
{
  for (auto pS : _pSystems) {
    pS->setVelToFluidVel(fVel);
  }
};

template<>
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::
setVelToFluidVel<descriptors::D3Q19<>>(
                                      SuperLatticeInterpPhysVelocity3D<double, descriptors::D3Q19<>>& fVel)
{
  for (auto pS : _pSystems) {
    pS->setVelToFluidVel(fVel);
  }
};

#ifndef OLB_PRECOMPILED
template<>
template<>
void SuperParticleSystem3D<double, Particle3D>::
setVelToFluidVel<descriptors::D3Q19<descriptors::FORCE>>(
      SuperLatticeInterpPhysVelocity3D<double, descriptors::D3Q19<descriptors::FORCE>>& fVel)
{
  for (auto pS : _pSystems) {
    pS->setVelToFluidVel(fVel);
  }
};

template<>
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::
setVelToFluidVel<descriptors::D3Q19<descriptors::FORCE>>(
      SuperLatticeInterpPhysVelocity3D<double, descriptors::D3Q19<descriptors::FORCE>>& fVel)
{
  for (auto pS : _pSystems) {
    pS->setVelToFluidVel(fVel);
  }
};
#endif


template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::setMagneticParticles(std::vector<double> dMoment,
                                                                             std::vector<double> vel, std::vector<double> aVel, std::vector<double> torque, double magnetisation)
{
    int i = 0;
    for (auto pS : _pSystems) {
        std::deque<MagneticParticle3D<double>*> particles = pS->getParticlesPointer();

        for (auto p : particles) {

            p->setMoment(dMoment);
            p->setVel(vel);
            p->setAVel(aVel);
            p->setTorque(torque);
            p->setMagnetisation(magnetisation);
            p->setID(i);
            i++;

        }
    }
}

template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::setMagneticParticles(std::vector<double> dMoment,
                                                                             std::vector<double> vel, std::vector<double> aVel, std::vector<double> torque, double magnetisation, int sActivity)
{
    int i = 0;
    for (auto pS : _pSystems) {
        std::deque<MagneticParticle3D<double>*> particles = pS->getParticlesPointer();

        for (auto p : particles) {

            p->setMoment(dMoment);
            p->setVel(vel);
            p->setAVel(aVel);
            p->setTorque(torque);
            p->setMagnetisation(magnetisation);
            p->setID(i);
            i++;
            p->setSActivity(sActivity);
        }
    }
}

template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::prepareAgglomerates()
{
    for (auto pS : _pSystems) {
        std::list<MagneticParticle3D<double>*> particlesList;
        pS->_Agglomerates.push_back(particlesList) ;
    }
}

template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::initAggloParticles()
{
    for (auto pS : _pSystems) {
        pS->initAggloParticles() ;
    }
}

template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::findAgglomerates(int iT, int itVtkOutputMagParticles)
{
    int pSi = 0;
    for (auto pS : _pSystems) {
        pS->findAgglomerates() ;

        if (iT % itVtkOutputMagParticles == 0) {

            clout << "Particlesystem number: " << pSi << std::endl;
            clout << "Number of non agglomerated particles" << ": " << pS->_Agglomerates[0].size() << std::endl;
            clout << "Number of agglomerated particles" << ": " << pS->size() - pS->_Agglomerates[0].size() << std::endl;
            clout << "Proportion of agglomeratet particles" << ": "
                  << double(pS->size() - pS->_Agglomerates[0].size()) / double(pS->size()) * 100. << "%" << std::endl;
            clout << "Number of agglomerates" << ": " << pS->_Agglomerates.size() - 1 << std::endl;
        }
        pSi++;
    }
}

// multiple collision models
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::simulate(double dT, std::set<int> sActivityOfParticle, bool scale)
{
    for (auto pS : _pSystems) {
        time_t delta = clock();
        if (pS->getIGeometry() == singleton::mpi().getRank()) {
            pS->_contactDetection->sort();
        }
        _stopSorting += clock() - delta;
        if (pS->getIGeometry() == singleton::mpi().getRank()) {
            pS->simulate(dT, sActivityOfParticle, scale);
            pS->computeBoundary();
        }
    }
    updateParticleDistribution();
}

template<>
bool SuperParticleSystem3D<double, MagneticParticle3D>::particleSActivityTest(int sActivity)
{
    for (auto pS : _pSystems) {
        for (auto p : pS->_particles) {
            if (p.getSActivity() == sActivity) {
                return false;
            }
        }
    }
    return true;
}

template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::setMagneticParticlesdMomRandom()
{

    for (auto pS : _pSystems) {
        std::deque<MagneticParticle3D<double>*> particles = pS->getParticlesPointer();

        for (auto p : particles) {
            std::vector<double> dMoment = { 0., 0., 0. };
            for (int i = 0; i < 3; i++) {
                dMoment[i] = rand() % (9 - (-9) + 1) + (-9);
            }

            double dMoment_norm = sqrt(pow(dMoment[0], 2.) + pow(dMoment[1], 2.) + pow(dMoment[2], 2.)) ;

            for (int i = 0; i < 3; i++) {
                dMoment[i] /= dMoment_norm ;
            }

            p->setMoment(dMoment);
        }
    }
}

template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::addParticle(
        IndicatorF3D<double>& ind, double mas, double rad, int no, int id,
        std::vector<double> vel, std::vector<double> dMoment, std::vector<double> aVel, std::vector<double> torque, double magnetisation,
        int sActivity)
{
    std::vector<double> pos(3, 0.);
    bool indic[1] = { false };

    no += globalNumOfParticles();
    while (globalNumOfParticles() < no) {
        pos[0] = ind.getMin()[0]
                 + (double) (rand() % 100000) / 100000. * (ind.getMax()[0] - ind.getMin()[0]);
        pos[1] = ind.getMin()[1]
                 + (double) (rand() % 100000) / 100000. * (ind.getMax()[1] - ind.getMin()[1]);
        pos[2] = ind.getMin()[2]
                 + (double) (rand() % 100000) / 100000. * (ind.getMax()[2] - ind.getMin()[2]);

#ifdef PARALLEL_MODE_MPI
        singleton::mpi().bCast(&*pos.begin(), 3);
#endif

        int x0, y0, z0, C;
        std::vector<int> locLat(4, 0);
        if (this->_cuboidGeometry.getFloorLatticeR(pos, locLat)) {
            C = locLat[0];
            if (this->_loadBalancer.rank(C) == singleton::mpi().getRank()) {
                x0 = locLat[1];
                y0 = locLat[2];
                z0 = locLat[3];
                if (_superGeometry.get(C, x0, y0, z0) == 1
                    && _superGeometry.get(C, x0, y0 + 1, z0) == 1
                    && _superGeometry.get(C, x0, y0, z0 + 1) == 1
                    && _superGeometry.get(C, x0, y0 + 1, z0 + 1) == 1
                    && _superGeometry.get(C, x0 + 1, y0, z0) == 1
                    && _superGeometry.get(C, x0 + 1, y0 + 1, z0) == 1
                    && _superGeometry.get(C, x0 + 1, y0, z0 + 1) == 1
                    && _superGeometry.get(C, x0 + 1, y0 + 1, z0 + 1) == 1
                    && ind(indic, &pos[0])) {
                    MagneticParticle3D<double> p(pos, vel, mas, rad, id, dMoment, aVel, torque, magnetisation, sActivity);
                    id++;
                    addParticle(p);
                }
            }
        }
    }
}

template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::addParticle(IndicatorF3D<double>& ind,  std::set<int>  material, double mas, double rad, int no, int id,
                                                                    std::vector<double> vel, std::vector<double> dMoment, std::vector<double> aVel, std::vector<double> torque, double magnetisation,
                                                                    int sActivity)

{
    std::vector<double> pos(3, 0.);
    bool indic[1] = { false };

    no += globalNumOfParticles();
    while (globalNumOfParticles() < no) {
        pos[0] = ind.getMin()[0]
                 + (double) (rand() % 100000) / 100000. * (ind.getMax()[0] - ind.getMin()[0]);
        pos[1] = ind.getMin()[1]
                 + (double) (rand() % 100000) / 100000. * (ind.getMax()[1] - ind.getMin()[1]);
        pos[2] = ind.getMin()[2]
                 + (double) (rand() % 100000) / 100000. * (ind.getMax()[2] - ind.getMin()[2]);

#ifdef PARALLEL_MODE_MPI
        singleton::mpi().bCast(&*pos.begin(), 3);
#endif

        int x0, y0, z0;
        std::vector<int> locLat(4, 0);
        if (this->_cuboidGeometry.getFloorLatticeR(pos, locLat)) {
            if (this->_loadBalancer.rank(locLat[0]) == singleton::mpi().getRank()) {
                x0 = locLat[1];
                y0 = locLat[2];
                z0 = locLat[3];
                if (_superGeometry.get(locLat[0], x0, y0, z0) == 1
                    && _superGeometry.get(locLat[0], x0, y0 + 1, z0) == 1
                    && _superGeometry.get(locLat[0], x0, y0, z0 + 1) == 1
                    && _superGeometry.get(locLat[0], x0, y0 + 1, z0 + 1) == 1
                    && _superGeometry.get(locLat[0], x0 + 1, y0, z0) == 1
                    && _superGeometry.get(locLat[0], x0 + 1, y0 + 1, z0) == 1
                    && _superGeometry.get(locLat[0], x0 + 1, y0, z0 + 1) == 1
                    && _superGeometry.get(locLat[0], x0 + 1, y0 + 1, z0 + 1) == 1
                    && ind(indic, &pos[0])) {
                    if (material.find(
                            _superGeometry.get(locLat[0], locLat[1], locLat[2], locLat[3]))
                        != material.end()
                        && material.find(
                            _superGeometry.get(locLat[0], locLat[1], locLat[2] + 1,
                                               locLat[3])) != material.end()
                        && material.find(
                            _superGeometry.get(locLat[0], locLat[1], locLat[2],
                                               locLat[3] + 1)) != material.end()
                        && material.find(
                            _superGeometry.get(locLat[0], locLat[1], locLat[2] + 1,
                                               locLat[3] + 1)) != material.end()
                        && material.find(
                            _superGeometry.get(locLat[0], locLat[1] + 1, locLat[2],
                                               locLat[3])) != material.end()
                        && material.find(
                            _superGeometry.get(locLat[0], locLat[1] + 1, locLat[2] + 1,
                                               locLat[3])) != material.end()
                        && material.find(
                            _superGeometry.get(locLat[0], locLat[1] + 1, locLat[2],
                                               locLat[3] + 1)) != material.end()
                        && material.find(
                            _superGeometry.get(locLat[0], locLat[1] + 1, locLat[2] + 1,
                                               locLat[3] + 1)) != material.end()) {

                        MagneticParticle3D<double> p(pos, vel, mas, rad, id, dMoment, aVel, torque, magnetisation, sActivity);
                        id++;
                        addParticle(p);
                    }
                }
            }
        }
    }
}

/*
// Magnetic particle type
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::addParticle(IndicatorF3D<double>& ind, double mas,
                                                                    double rad, int no, int id,
                                                                    std::vector<double> vel, std::vector<double> dMoment, std::vector<double> aVel,
                                                                    std::vector<double> torque, double magnetisation, int sActivity);
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::addParticle(IndicatorF3D<double>& ind, double mas,
                                                                    double rad, int no, int id,
                                                                    std::vector<double> vel, std::vector<double> dMoment, std::vector<double> aVel,
                                                                    std::vector<double> torque, double magnetisation, int sActivity);
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::addParticle(IndicatorF3D<double>& ind,
                                                                    std::set<int>  material, double mas, double rad, int no, int id,
                                                                    std::vector<double> vel, std::vector<double> dMoment, std::vector<double> aVel,
                                                                    std::vector<double> torque, double magnetisation, int sActivity);

template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::setMagneticParticlesdMomRandom();
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::setMagneticParticles(std::vector<double> dMoment,
                                                                             std::vector<double> vel, std::vector<double> aVel,
                                                                             std::vector<double> torque, double magnetisation, int sActivity);
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::prepareAgglomerates();
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::initAggloParticles();
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::findAgglomerates(int iT, int itVtkOutputMagParticles);
template<>
bool SuperParticleSystem3D<double, MagneticParticle3D>::particleSActivityTest(int sActivity);
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::simulate(double dT, std::set<int> sActivityOfFreeParticle, bool scale) ;
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::simulateWithTwoWayCoupling_Mathias ( double dT,
                                                                                             ForwardCouplingModel<double,MagneticParticle3D>& forwardCoupling,
                                                                                             BackCouplingModel<double,MagneticParticle3D>& backCoupling,
                                                                                             int material, int subSteps, bool resetExternalField, bool scale );
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::simulateWithTwoWayCoupling_Davide ( double dT,
                                                                                            ForwardCouplingModel<double,MagneticParticle3D>& forwardCoupling,
                                                                                            BackCouplingModel<double,MagneticParticle3D>& backCoupling,
                                                                                            int material, int subSteps, bool resetExternalField, bool scale );
*/
}  // namespace olb
