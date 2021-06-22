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

#include "particle3D.h"
#include "particle3D.hh"
#include "particleSpecializations/particleSpecializations3D.h"
#include "particleSpecializations/particleSpecializations3D.hh"
#include "particleSystem3D.h"
#include "particleSystem3D.hh"
#include "dynamics/latticeDescriptors.h"

#include "functors/lattice/latticeInterpPhysVelocity3D.h"

namespace olb {


template class SimulateParticles<double,Particle3D>;
template class ParticleSystem3D<double,Particle3D>;

template class SimulateParticles<double, MagneticParticle3D>;
template class ParticleSystem3D<double, MagneticParticle3D>;

template<>
template<>
void ParticleSystem3D<double,Particle3D>::
setVelToFluidVel<descriptors::D3Q19<>>(
                                      SuperLatticeInterpPhysVelocity3D<double, descriptors::D3Q19<>>& fVel)
{
  for (auto& p : _particles) {
    if (p.getActive()) {
      fVel(&p.getVel()[0], &p.getPos()[0], p.getCuboid());
    }
  }
};

template<>
template<>
void ParticleSystem3D<double,MagneticParticle3D>::
setVelToFluidVel<descriptors::D3Q19<>>(
                                      SuperLatticeInterpPhysVelocity3D<double, descriptors::D3Q19<>>& fVel)
{
  for (auto& p : _particles) {
    if (p.getActive()) {
      fVel(&p.getVel()[0], &p.getPos()[0], p.getCuboid());
    }
  }
};

#ifndef OLB_PRECOMPILED
template<>
template<>
void ParticleSystem3D<double,Particle3D>::
setVelToFluidVel<descriptors::D3Q19<descriptors::FORCE>>(
      SuperLatticeInterpPhysVelocity3D<double, descriptors::D3Q19<descriptors::FORCE>>& fVel)
{
  for (auto& p : _particles) {
    if (p.getActive()) {
      fVel(&p.getVel()[0], &p.getPos()[0], p.getCuboid());
    }
  }
};

template<>
template<>
void ParticleSystem3D<double,MagneticParticle3D>::
setVelToFluidVel<descriptors::D3Q19<descriptors::FORCE>>(
      SuperLatticeInterpPhysVelocity3D<double, descriptors::D3Q19<descriptors::FORCE>>& fVel)
{
  for (auto& p : _particles) {
    if (p.getActive()) {
      fVel(&p.getVel()[0], &p.getPos()[0], p.getCuboid());
    }
  }
};
#endif

template<>
void ParticleSystem3D<double, MagneticParticle3D>::explicitEuler(double dT, bool scale)
{
    double maxDeltaR = _superGeometry.getCuboidGeometry().getMaxDeltaR();
    double maxFactor = double();

    for (auto& p : _particles) {

        if (p.getActive()) {

            if (p.getSActivity() == 3) {continue;}

            for (int i = 0; i < 3; i++) {
                p.getVel()[i] += p.getForce()[i] * p.getInvAddedMass() * dT;
                p.getPos()[i] += p.getVel()[i] * dT;

                // computation of direction depending maxFactor to scale velocity value
                // to keep velocity small enough for simulation
                if (fabs(p.getVel()[i]) > fabs(maxDeltaR / dT)) {
                    maxFactor = std::max(maxFactor, fabs(p.getVel()[i] / maxDeltaR * dT));
                }
            }
            // scaling of velocity values
            // if particles are too fast, e.g. material boundary can not work anymore
            if ( !util::nearZero(maxFactor) && scale) {
                std::cout << "particle velocity is scaled because of reached limit"
                          << std::endl;
                for (int i = 0; i < 3; i++) {
                    p.getPos()[i] -= p.getVel()[i] * dT; // position set back
                    p.getVel()[i] /= maxFactor; // scale velocity value
                    p.getPos()[i] += p.getVel()[i] * dT;
                }
            }

            // Change sActivity of particle in dependence of its position in the geometry
            // if ((p.getPos()[0] > 0.00015) && (p.getSActivity() == 0)) {
            //   p.setSActivity(2);
            // }

            // }

            // if particles are too fast, e.g. material boundary can not work anymore
//#ifdef OLB_DEBUG
//      if (p.getVel()[i] * dT > _superGeometry.getCuboidGeometry().getMaxDeltaR()) {
//        std::cout << " PROBLEM: particle speed too high rel. to delta of "
//                  "lattice: "<< std::endl;
//        std::cout << "p.getVel()[i]*dT: " << i <<" "<< p.getVel()[i] * dT;
//        std::cout << "MaxDeltaR(): " <<
//                  _superGeometry.getCuboidGeometry().getMaxDeltaR() << std::endl;
//        exit(-1);
//      }
//#endif

        }
    }
}

// multiple collision models
template<>
void ParticleSystem3D<double, MagneticParticle3D>::explicitEuler(double dT, std::set<int> sActivityOfParticle, bool scale)
{
    double maxDeltaR = _superGeometry.getCuboidGeometry().getMaxDeltaR();
    double maxFactor = double();

    for (auto& p : _particles) {

        if (p.getActive()) {

            if (p.getSActivity() == 3) {continue;}

            bool b = false;
            for (auto sA : sActivityOfParticle) {
                if (p.getSActivity() == sA) {b = true; break;}
            }
            if (b == false) {continue;}

            for (int i = 0; i < 3; i++) {
                p.getVel()[i] += p.getForce()[i] * p.getInvAddedMass() * dT;
                p.getPos()[i] += p.getVel()[i] * dT;

                // computation of direction depending maxFactor to scale velocity value
                // to keep velocity small enough for simulation
                if (fabs(p.getVel()[i]) > fabs(maxDeltaR / dT)) {
                    maxFactor = std::max(maxFactor, fabs(p.getVel()[i] / maxDeltaR * dT));
                }
            }
            // scaling of velocity values
            // if particles are too fast, e.g. material boundary can not work anymore
            if ( !util::nearZero(maxFactor) && scale) {
                std::cout << "particle velocity is scaled because of reached limit"
                          << std::endl;
                for (int i = 0; i < 3; i++) {
                    p.getPos()[i] -= p.getVel()[i] * dT; // position set back
                    p.getVel()[i] /= maxFactor; // scale velocity value
                    p.getPos()[i] += p.getVel()[i] * dT;
                }
            }

            // Change sActivity of particle in dependence of its position in the geometry
            // if ((p.getPos()[0] > 0.00015) && (p.getSActivity() == 0)) {
            //   p.setSActivity(2);
            // }

            // }

            // if particles are too fast, e.g. material boundary can not work anymore
//#ifdef OLB_DEBUG
//      if (p.getVel()[i] * dT > _superGeometry.getCuboidGeometry().getMaxDeltaR()) {
//        std::cout << " PROBLEM: particle speed too high rel. to delta of "
//                  "lattice: "<< std::endl;
//        std::cout << "p.getVel()[i]*dT: " << i <<" "<< p.getVel()[i] * dT;
//        std::cout << "MaxDeltaR(): " <<
//                  _superGeometry.getCuboidGeometry().getMaxDeltaR() << std::endl;
//        exit(-1);
//      }
//#endif

        }
    }
}

//template<typename T, template<typename U> class PARTICLETYPE>
//void ParticleSystem3D<T, PARTICLETYPE>::integrateTorque(T dT)
//{
//  for (auto& p : _particles) {
//    if (p.getActive()) {
//      for (int i = 0; i < 3; i++) {
//        p.getAVel()[i] += p.getTorque()[i] * 1.
//                          / (2. / 5. * p.getMass() * std::pow(p.getRad(), 2)) * dT;
//      }
//    }
//  }
//}

//template<typename T, template<typename U> class PARTICLETYPE>
//void ParticleSystem3D<T, PARTICLETYPE>::integrateTorqueMag(T dT) {
////template<typename T>
////void ParticleSystem3D<T, MagneticParticle3D>::integrateTorqueMag(T dT) {
//  for (auto& p : _particles) {
//    if (p.getActive()) {
//      Vector<T, 3> deltaAngle;
//      T angle;
//      T epsilon = std::numeric_limits<T>::epsilon();
//      for (int i = 0; i < 3; i++) {
//        // change orientation due to torque moments
//        deltaAngle[i] = (5. * p.getTorque()[i] * dT * dT) / (2. * p.getMass() * std::pow(p.getRad(), 2));
//        // apply change in angle to dMoment vector
//      }
//      angle = norm(deltaAngle);
//      if (angle > epsilon) {
//        //Vector<T, 3> axis(deltaAngle);
//        //  Vector<T, 3> axis(T(), T(), T(1));
//        //axis.normalize();
//        std::vector<T> null(3, T());
//
//        //RotationRoundAxis3D<T, S> rotRAxis(p.getPos(), fromVector3(axis), angle);
//        RotationRoundAxis3D<T, S> rotRAxis(null, fromVector3(deltaAngle), angle);
//        T input[3] = {p.getMoment()[0], p.getMoment()[1], p.getMoment()[2]};
//        Vector<T, 3> in(input);
//        T output[3] = {T(), T(), T()};
//        rotRAxis(output, input);
//        Vector<T, 3> out(output);
////        std::vector<T> mainRefVec(3, T());
////        mainRefVec[0] = 1.;
////        std::vector<T> secondaryRefVec(3, T());
////        secondaryRefVec[2] = 1.;
////        AngleBetweenVectors3D<T, T> checkAngle(mainRefVec, secondaryRefVec);
////        T angle[1];
////        checkAngle(angle, input);
//        std::cout<< "|moment|_in: " << in.norm() << ", |moment|_out: " << out.norm()
//                 << ", | |in| - |out| |: " << fabs(in.norm() - out.norm())
//                 /*<< " Angle: " << angle[0] */<< std::endl;
//        p.getMoment()[0] = output[0];
//        p.getMoment()[1] = output[1];
//        p.getMoment()[2] = output[2];
//      }
//    }
//  }
//}

//template<typename T, template<typename U> class PARTICLETYPE>
//void ParticleSystem3D<T, PARTICLETYPE>::implicitEuler(T dT, AnalyticalF<3,T,T>& getvel) {
//  _activeParticles = 0;
//  for (auto& p : _particles) {
//    if(p.getActive()) {
//      std::vector<T> fVel = getvel(p._pos);
////      std::vector<T> vel = p.getVel();
//      std::vector<T> pos = p.getPos();
//      T C = 6.* M_PI * p._rad * this->_converter.getCharNu() * this->_converter.getCharRho()* dT/p.getMass();
//      for (int i = 0; i<3; i++) {
//        p._vel[i] = (p._vel[i]+C*fVel[i]) / (1+C);
//        p._pos[i] += p._vel[i] * dT;
//      }
////      p.setVel(vel);
////      p.setPos(pos);
//      checkActive(p);
////      cout << "C: " << C << std::endl;
//    }
////    cout << "pos: " << p._pos[0] << " " << p._pos[1] << " " << p._pos[2] << " " << p._vel[0] << " " << p._vel[1] << " " << p._vel[2]<< std::endl; //"\n force: " << p._force[0] << " " << p._force[1] << " " << p._force[2]<< "\n " << std::endl;
//  }
//}

/*
template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::predictorCorrector1(T dT)
{
  std::vector<T> vel;
  std::vector<T> pos;
  std::vector<T> frc;
  for (auto& p : _particles) {
    if (p.getActive()) {
      vel = p.getVel();
      p.setVel(vel, 1);
      pos = p.getPos();
      p.setVel(pos, 2);
      frc = p.getForce();
      p.setForce(frc, 1);
      for (int i = 0; i < 3; i++) {
        vel[i] += p._force[i] / p.getMass() * dT;
        pos[i] += vel[i] * dT;
      }
      p.setVel(vel);
      p.setPos(pos);
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::predictorCorrector2(T dT)
{
  std::vector<T> vel;
  std::vector<T> pos;
  std::vector<T> frc;
  for (auto& p : _particles) {
    if (p.getActive()) {
      vel = p.getVel(1);
      pos = p.getVel(2);
      for (int i = 0; i < 3; i++) {
        vel[i] += dT * .5 * (p.getForce()[i] + p.getForce(1)[i]) / p.getMass();
        pos[i] += vel[i] * dT;
      }
      p.setVel(vel);
      p.setPos(pos);
    }
  }
}
*/

/*
template<typename T, template<typename U> class PARTICLETYPE>
void ParticleSystem3D<T, PARTICLETYPE>::adamBashforth4(T dT)
{
  for (auto& p : _particles) {
    if (p.getActive()) {
      std::vector<T> vel = p.getVel();
      std::vector<T> pos = p.getPos();
      for (int i = 0; i < 3; i++) {
        vel[i] += dT / p.getMas()
                  * (55. / 24. * p.getForce()[i] - 59. / 24. * p.getForce(1)[i]
                     + 37. / 24. * p.getForce(2)[i] - 9. / 24. * p.getForce(3)[i]);
      }
      p.rotAndSetVel(vel);
      for (int i = 0; i < 3; i++) {
        pos[i] += dT
                  * (55. / 24. * p.getVel()[i] - 59. / 24. * p.getVel(1)[i]
                     + 37. / 24. * p.getVel(2)[i] - 3. / 8. * p.getVel(3)[i]);
      }
      p.setPos(pos);
    }
  }
}
*/

template<>
void ParticleSystem3D<double, MagneticParticle3D>::resetMag()
{
    typename std::deque<MagneticParticle3D<double> >::iterator p;
    int pInt = 0;
    for (p = _particles.begin(); p != _particles.end(); ++p, ++pInt) {

        if (p->getActive()) {

            p->resetForce();
            p->resetTorque();
        }
    }
}

// multiple collision models
template<>
void ParticleSystem3D<double, MagneticParticle3D>::resetMag(std::set<int> sActivityOfParticle)
{
    typename std::deque<MagneticParticle3D<double> >::iterator p;
    int pInt = 0;
    for (p = _particles.begin(); p != _particles.end(); ++p, ++pInt) {

        if (p->getSActivity() == 3) {continue;}

        if (p->getActive()) {
            bool b = false;
            for (auto sA : sActivityOfParticle) {
                if (p->getSActivity() == sA) {b = true; break;}
            }
            if (b == false) {continue;}
            p->resetForce();
            p->resetTorque();
        }
    }
}


template<>
void ParticleSystem3D<double, MagneticParticle3D>::computeForce()
{
    typename std::deque<MagneticParticle3D<double> >::iterator p;
    int pInt = 0;
    for (p = _particles.begin(); p != _particles.end(); ++p, ++pInt) {
        if (p->getActive()) {

            for (auto f : _forces) {
                f->applyForce(p, pInt, *this);
            }
        }
    }
}

// multiple collision models
template<>
void ParticleSystem3D<double, MagneticParticle3D>::computeForce(std::set<int> sActivityOfParticle)
{
    typename std::deque<MagneticParticle3D<double> >::iterator p;
    int pInt = 0;
    for (p = _particles.begin(); p != _particles.end(); ++p, ++pInt) {

        if (p->getActive()) {

            if (p->getSActivity() == 3) {continue;}

            bool b = false;
            for (auto sA : sActivityOfParticle) {
                if (p->getSActivity() == sA) {b = true; break;}
            }
            if (b == false) {continue;}

            for (auto f : _forces) {
                // f->applyForce(p, p->getID(), *this);
                f->applyForce(p, pInt, *this);
            }
        }
    }
}

/* Original MagDM damping intern
template<>
void ParticleSystem3D<double, MagneticParticle3D>::integrateTorqueMag(double dT)
{
  for (auto& p : _particles) {
    Vector<double, 3> deltaAngle;
    double angle;
    double epsilon = std::numeric_limits<double>::epsilon();
    double damping = std::pow((1. - p.getADamping()), dT);
    for (int i = 0; i < 3; i++) {
      p.getAVel()[i] += (5. * (p.getTorque()[i]) * dT) / (2.  * p.getMass() * std::pow(p.getRad(), 2));
      p.getAVel()[i] *= damping;
      deltaAngle[i] = p.getAVel()[i] * dT;

    angle = norm(deltaAngle);
    if (angle > epsilon) {
      std::vector<double> null(3, double());

      RotationRoundAxis3D<double, double> rotRAxis(null, util::fromVector3(deltaAngle), angle);
      double input[3] = {p.getMoment()[0], p.getMoment()[1], p.getMoment()[2]};
      Vector<double, 3> in(input);
      double output[3] = {double(), double(), double()};
      rotRAxis(output, input);
      Vector<double, 3> out(output);
      // renormalize output
      if (out.norm() > epsilon) {
        out = (1. / out.norm()) * out;
      }

      p.getMoment()[0] = out[0];
      p.getMoment()[1] = out[1];
      p.getMoment()[2] = out[2];
    }
  }
}
}
*/

template<>
void ParticleSystem3D<double, MagneticParticle3D>::integrateTorqueMag(double dT)
{
    for (auto& p : _particles) {

        Vector<double, 3> deltaAngle;
        double angle;
        double epsilon = std::numeric_limits<double>::epsilon();
        for (int i = 0; i < 3; i++) {
            p.getAVel()[i] += (5. * (p.getTorque()[i]) * dT) / (2.  * p.getMass() * std::pow(p.getRad(), 2));
            deltaAngle[i] = p.getAVel()[i] * dT;
            angle = norm(deltaAngle);

            if (angle > epsilon) {
                std::vector<double> null(3, double());
                RotationRoundAxis3D<double, double> rotRAxis(null, util::fromVector3(deltaAngle), angle);
                double input[3] = {p.getMoment()[0], p.getMoment()[1], p.getMoment()[2]};
                Vector<double, 3> in(input);
                double output[3] = {double(), double(), double()};
                rotRAxis(output, input);
                Vector<double, 3> out(output);
                // renormalize output
                if (norm(out) > epsilon) {
                    out = normalize(out);
                }

                p.getMoment()[0] = out[0];
                p.getMoment()[1] = out[1];
                p.getMoment()[2] = out[2];
            }
        }
    }
}

// multiple collision models
template<>
void ParticleSystem3D<double, MagneticParticle3D>::integrateTorqueMag(double dT, std::set<int> sActivityOfParticle)
{

    for (auto& p : _particles) {

        if (p.getSActivity() == 3) {continue;}

        bool b = false;
        for (auto sA : sActivityOfParticle) {
            if (p.getSActivity() == sA) {b = true; break;}
        }
        if (b == false) {continue;}

        Vector<double, 3> deltaAngle;
        double angle;
        double epsilon = std::numeric_limits<double>::epsilon();
        for (int i = 0; i < 3; i++) {
            p.getAVel()[i] += (5. * (p.getTorque()[i]) * dT) / (2.  * p.getMass() * std::pow(p.getRad(), 2));
            deltaAngle[i] = p.getAVel()[i] * dT;
            angle = norm(deltaAngle);

            if (angle > epsilon) {
                std::vector<double> null(3, double());
                RotationRoundAxis3D<double, double> rotRAxis(null, util::fromVector3(deltaAngle), angle);
                double input[3] = {p.getMoment()[0], p.getMoment()[1], p.getMoment()[2]};
                Vector<double, 3> in(input);
                double output[3] = {double(), double(), double()};
                rotRAxis(output, input);
                Vector<double, 3> out(output);
                // renormalize output
                if (norm(out) > epsilon) {
                    out = normalize(out);
                }

                p.getMoment()[0] = out[0];
                p.getMoment()[1] = out[1];
                p.getMoment()[2] = out[2];
            }
        }
    }
}

template<>
void ParticleSystem3D<double, MagneticParticle3D>::setOverlapZero()
{

    typename std::deque<MagneticParticle3D<double> >::iterator p;
    int pInt = 0;
    for (p = _particles.begin(); p != _particles.end(); ++p, ++pInt) {

        std::vector<std::pair<size_t, double>> ret_matches;
        // kind of contactDetection has to be chosen in application
        getContactDetection()->getMatches(pInt, ret_matches);

        MagneticParticle3D<double>* p2 = NULL;

        // iterator walks through number of neighbored particles = ret_matches
        for (const auto& it : ret_matches) {

            if (!util::nearZero(it.second)) {
                p2 = &(_particles.at(it.first)); //p2 = &pSys[it.first];

                if ((p2->getRad() + p->getRad()) > (sqrt(it.second))) {

                    // overlap
                    double overlap = (p2->getRad() + p->getRad()) - sqrt(it.second);

                    //conVec: vector from particle1 to particle2
                    Vector<double, 3> conVec(0., 0., 0.) ;
                    for (int i = 0; i <= 2; i++) {

                        conVec[i] = p2->getPos()[i] - p->getPos()[i];
                    }
                    Vector<double, 3> conVecNormalized(conVec) ;
                    normalize(conVecNormalized) ;

                    double dpos[3] = {double(0), double(0), double(0) } ;

                    // Both particles are not deposited (sActivity = 3)
                    if ((p->getSActivity() != 3) && (p2->getSActivity() != 3)) {

                        for (int i = 0; i <= 2; i++) {

                            dpos[i] = conVecNormalized[i] * 0.5 * overlap ;
                            p->getPos()[i] -= 1.* dpos[i];
                            p2->getPos()[i] += 1.* dpos[i];
                        }
                        if ((p->getSActivity() == 2) || (p->getSActivity() == 2)) {
                            p->setSActivity(2);
                            p2->setSActivity(2);
                        }
                        continue;
                    }

                    // Particle 1 is deposited (sActivity = 3) and Particle 2 is not
                    if ((p->getSActivity() != 3) && (p2->getSActivity() == 3)) {

                        for (int i = 0; i <= 2; i++) {

                            dpos[i] = conVecNormalized[i] * 0.5 * overlap ;
                            p->getPos()[i] -= 2. * dpos[i];
                        }
                        p->setSActivity(2) ;
                    }

                    // Particle 2 is deposited (sActivity = 3) and Particle 1 is not
                    if ((p->getSActivity() == 3) && (p2->getSActivity() != 3)) {

                        for (int i = 0; i <= 2; i++) {

                            dpos[i] = conVecNormalized[i] * 0.5 * overlap ;
                            p2->getPos()[i] += 2. * dpos[i];
                        }
                        p2->setSActivity(2) ;
                    }

                    // Both particles are not deposited (sActivity = 3)
                    if ((p->getSActivity() == 3) && (p2->getSActivity() == 3)) {

                        for (int i = 0; i <= 2; i++) {

                            dpos[i] = conVecNormalized[i] * 0.5 * overlap ;
                            p->getPos()[i] -= dpos[i];
                            p2->getPos()[i] += dpos[i];
                        }
                    }
                }
            }
        }
    }
}

template<>
void ParticleSystem3D<double, MagneticParticle3D>::setOverlapZeroForCombinationWithMechContactForce()
{

    typename std::deque<MagneticParticle3D<double> >::iterator p;
    int pInt = 0;

    for (p = _particles.begin(); p != _particles.end(); ++p, ++pInt) {

        if (p->getSActivity() > 1) { continue; }

        std::vector<std::pair<size_t, double>> ret_matches;
        // kind of contactDetection has to be chosen in application
        getContactDetection()->getMatches(pInt, ret_matches);

        MagneticParticle3D<double>* p2 = NULL;
        // iterator walks through number of neighbored particles = ret_matches
        for (const auto& it : ret_matches) {
            if (!util::nearZero(it.second)) {
                p2 = &(_particles.at(it.first)); //p2 = &pSys[it.first];

                if ((p2->getRad() + p->getRad()) > (sqrt(it.second))) {
                    // overlap
                    double overlap = (p2->getRad() + p->getRad()) - sqrt(it.second);

                    //conVec: vector from particle1 to particle2
                    Vector<double, 3> conVec(0., 0., 0.) ;
                    for (int i = 0; i <= 2; i++) {

                        conVec[i] = p2->getPos()[i] - p->getPos()[i];
                    }
                    Vector<double, 3> conVecNormalized(conVec) ;
                    normalize(conVecNormalized) ;

                    double dpos[3] = {double(0), double(0), double(0) } ;

                    // Both particles are not deposited
                    if (overlap > p2->getRad() + p->getRad()) {

                        if (p2->getSActivity() == 1)  {

                            for (int i = 0; i <= 2; i++) {

                                dpos[i] = conVecNormalized[i] * 0.5 * overlap ;
                                p->getPos()[i] -= dpos[i] * 1. ;
                                p2->getPos()[i] += dpos[i] * 1. ;
                            }
                            continue;
                        }
                    }

                        // Particle 2 is out of the sphere of influence of setOverlapZeroForCombinationWithMechContactForce()
                        // Particle 1 is transferred to the influence of the mechanic contact force
                    else {
                        for (int i = 0; i <= 2; i++) {

                            dpos[i] = conVecNormalized[i] * 0.5 * overlap ;
                            p->getPos()[i] -= 2 * dpos[i];
                        }
                        p->setSActivity(2);
                        continue;
                    }
                }
            }
        }
    }
}

template<>
void ParticleSystem3D<double, MagneticParticle3D>::partialElasticImpact(double restitutionCoeff)
{
    typename std::deque<MagneticParticle3D<double> >::iterator p;
    int pInt = 0;

    for (p = _particles.begin(); p != _particles.end(); ++p, ++pInt) {

        std::vector<std::pair<size_t, double>> ret_matches;
        // kind of contactDetection has to be chosen in application
        getContactDetection()->getMatches(pInt, ret_matches);

        MagneticParticle3D<double>* p2 = NULL;

        // iterator walks through number of neighbored particles = ret_matches
        for (const auto& it : ret_matches) {

            if (!util::nearZero(it.second)) {
                p2 = &(_particles.at(it.first));

                if ((p2->getRad() + p->getRad()) > (sqrt(it.second))) {

                    // overlap
                    double overlap = (p2->getRad() + p->getRad()) - sqrt(it.second);

                    //conVec: vector from particle1 to particle2
                    Vector<double, 3> conVec(0., 0., 0.) ;
                    for (int i = 0; i <= 2; i++) {

                        conVec[i] = p2->getPos()[i] - p->getPos()[i];
                    }
                    Vector<double, 3> conVecNormalized(conVec) ;
                    normalize(conVecNormalized) ;

                    // Particle velocities before collision
                    Vector<double, 3> vel1bc = {p->getVel()} ;
                    Vector<double, 3> vel2bc = {p2->getVel()} ;

                    // Normal and tangential particle velocities before collision
                    Vector<double, 3> velN1(0., 0., 0.) ;
                    Vector<double, 3> velT1(0., 0., 0.) ;
                    Vector<double, 3> velN2(0., 0., 0.) ;
                    Vector<double, 3> velT2(0., 0., 0.) ;

                    // Particle velocities after collision
                    Vector<double, 3> vel1ac(0., 0., 0.) ;
                    Vector<double, 3> vel2ac(0., 0., 0.) ;

                    // Restitution coeffizient
                    double Cr = restitutionCoeff;

                    velN1 = (conVecNormalized * vel1bc) * conVecNormalized;
                    velN2 = (conVecNormalized * vel2bc) * conVecNormalized;
                    velT1 = vel1bc - velN1;
                    velT2 = vel2bc - velN2;
                    vel1ac = (Cr * p2->getMass() * ( velN2 - velN1) + (p2->getMass() * velN2) + (p->getMass() * velN1)) * (1. / (p->getMass() + p2->getMass())) ;
                    vel2ac = (Cr * p->getMass() * ( velN1 - velN2) + (p2->getMass() * velN2) + (p->getMass() * velN1)) * (1. / (p->getMass() + p2->getMass())) ;

                    double dpos[3] = {double(0), double(0), double(0) } ;

                    // Both particles are not deposited (sActivity = 3)
                    if ((p->getSActivity() != 3) && (p2->getSActivity() != 3)) {

                        for (int i = 0; i <= 2; i++) {

                            dpos[i] = conVecNormalized[i] * 0.5 * overlap ;
                            p->getPos()[i] -= dpos[i] * 1. ;
                            p->getVel()[i] = vel1ac[i] + velT1[i] ;
                            p2->getPos()[i] += dpos[i] * 1. ;
                            p2->getVel()[i] = vel2ac[i] + velT2[i] ;
                        }
                        if ((p->getSActivity() == 2) || (p->getSActivity() == 2)) {
                            p->setSActivity(2);
                            p2->setSActivity(2);
                        }
                        continue;
                    }

                    // Particle 1 is deposited (sActivity = 3) and Particle 2 is not
                    if ((p->getSActivity() != 3) && (p2->getSActivity() == 3)) {

                        for (int i = 0; i <= 2; i++) {

                            dpos[i] = conVecNormalized[i] * 0.5 * overlap ;
                            p->getPos()[i] -= 2. * dpos[i];
                            p->getVel()[i] = -1. * p->getVel()[i];
                        }
                        p->setSActivity(2) ;
                        continue;
                    }

                    // Particle 2 is deposited (sActivity = 3) and Particle 1 is not
                    if ((p->getSActivity() == 3) && (p2->getSActivity() != 3)) {

                        for (int i = 0; i <= 2; i++) {

                            dpos[i] = conVecNormalized[i] * 0.5 * overlap ;
                            p2->getPos()[i] += 2. * dpos[i];
                            p2->getVel()[i] = -1. * p2->getVel()[i];
                        }
                        p2->setSActivity(2) ;
                        continue;
                    }

                    // Both particles are deposited (sActivity = 3)
                    if ((p->getSActivity() == 3) && (p2->getSActivity() == 3)) {

                        for (int i = 0; i <= 2; i++) {

                            dpos[i] = conVecNormalized[i] * 0.5 * overlap ;
                            p->getPos()[i] -= dpos[i];
                            p2->getPos()[i] += dpos[i];
                        }
                        continue;
                    }
                }
            }
        }
    }
}

template<>
void ParticleSystem3D<double, MagneticParticle3D>::partialElasticImpactV2(double restitutionCoeff)
{
    typename std::deque<MagneticParticle3D<double> >::iterator p;
    int pInt = 0;

    for (p = _particles.begin(); p != _particles.end(); ++p, ++pInt) {

        std::vector<std::pair<size_t, double>> ret_matches;
        // kind of contactDetection has to be chosen in application
        getContactDetection()->getMatches(pInt, ret_matches);

        MagneticParticle3D<double>* p2 = NULL;

        // iterator walks through number of neighbored particles = ret_matches
        if (ret_matches.size() == 1) {continue;}
        double minDist = 0. ;
        int xPInt = 0 ;

        for (auto a : ret_matches) {
            if (util::nearZero(a.second)) {continue;}
            if (minDist == 0) {minDist = a.second; xPInt = a.first;}
            else {
                if (a.second < minDist) {minDist = a.second; xPInt = a.first;}
                continue;
            }
        }

        p2 = &(_particles.at(xPInt));

        if ((p2->getRad() + p->getRad()) > (sqrt(minDist))) {

            // overlap
            double overlap = (p2->getRad() + p->getRad()) - sqrt(minDist);

            //conVec: vector from particle 1 to particle 2
            Vector<double, 3> conVec(0., 0., 0.) ;
            for (int i = 0; i <= 2; i++) {

                conVec[i] = p2->getPos()[i] - p->getPos()[i];
            }
            Vector<double, 3> conVecNormalized(conVec) ;
            normalize(conVecNormalized) ;

            // Particle velocities before collision
            Vector<double, 3> vel1bc = {p->getVel()} ;
            Vector<double, 3> vel2bc = {p2->getVel()} ;

            // Normal and tangential particle velocities before collision
            Vector<double, 3> velN1(0., 0., 0.) ;
            Vector<double, 3> velT1(0., 0., 0.) ;
            Vector<double, 3> velN2(0., 0., 0.) ;
            Vector<double, 3> velT2(0., 0., 0.) ;

            // Particle velocities after collision
            Vector<double, 3> vel1ac(0., 0., 0.) ;
            Vector<double, 3> vel2ac(0., 0., 0.) ;

            // Restitution coeffizient
            double Cr = restitutionCoeff;

            velN1 = (conVecNormalized * vel1bc) * conVecNormalized;
            velN2 = (conVecNormalized * vel2bc) * conVecNormalized;
            velT1 = vel1bc - velN1;
            velT2 = vel2bc - velN2;
            vel1ac = (Cr * p2->getMass() * ( velN2 - velN1) + (p2->getMass() * velN2) + (p->getMass() * velN1)) * (1. / (p->getMass() + p2->getMass())) ;
            vel2ac = (Cr * p->getMass() * ( velN1 - velN2) + (p2->getMass() * velN2) + (p->getMass() * velN1)) * (1. / (p->getMass() + p2->getMass())) ;

            double dpos[3] = {double(0), double(0), double(0) } ;

            // Both particles are not deposited (sActivity = 3)
            if ((p->getSActivity() != 3) && (p2->getSActivity() != 3)) {

                for (int i = 0; i <= 2; i++) {

                    dpos[i] = conVecNormalized[i] * 0.5 * overlap ;
                    p->getPos()[i] -= dpos[i] * 1. ;
                    p->getVel()[i] = vel1ac[i] + velT1[i] ;
                    p2->getPos()[i] += dpos[i] * 1. ;
                    p2->getVel()[i] = vel2ac[i] + velT2[i] ;
                }
                if ((p->getSActivity() == 2) || (p->getSActivity() == 2)) {
                    p->setSActivity(2);
                    p2->setSActivity(2);
                }
                continue;
            }

            // Particle 1 is deposited (sActivity = 3) and Particle 2 is not
            if ((p->getSActivity() != 3) && (p2->getSActivity() == 3)) {

                for (int i = 0; i <= 2; i++) {

                    dpos[i] = conVecNormalized[i] * 0.5 * overlap ;
                    p->getPos()[i] -= 2. * dpos[i];
                    p->getVel()[i] = -1. * p->getVel()[i];
                }
                p->setSActivity(2) ;
                continue;
            }

            // Particle 2 is deposited (sActivity = 3) and Particle 1 is not
            if ((p->getSActivity() == 3) && (p2->getSActivity() != 3)) {

                for (int i = 0; i <= 2; i++) {

                    dpos[i] = conVecNormalized[i] * 0.5 * overlap ;
                    p2->getPos()[i] += 2. * dpos[i];
                    p2->getVel()[i] = -1. * p2->getVel()[i];
                }
                p2->setSActivity(2) ;
                continue;
            }

            // Both particles are not deposited (sActivity = 3)
            if ((p->getSActivity() == 3) && (p2->getSActivity() == 3)) {

                for (int i = 0; i <= 2; i++) {

                    dpos[i] = conVecNormalized[i] * 0.5 * overlap ;
                    p->getPos()[i] -= dpos[i];
                    p2->getPos()[i] += dpos[i];
                }
                continue;
            }
        }
    }
}

template<>
void ParticleSystem3D<double, MagneticParticle3D>::partialElasticImpactForCombinationWithMechContactForce(double restitutionCoeff)
{
    typename std::deque<MagneticParticle3D<double> >::iterator p;
    int pInt = 0;

    for (p = _particles.begin(); p != _particles.end(); ++p, ++pInt) {

        if (p->getSActivity() > 1) { continue; }

        std::vector<std::pair<size_t, double>> ret_matches;
        // kind of contactDetection has to be chosen in application
        getContactDetection()->getMatches(pInt, ret_matches);

        MagneticParticle3D<double>* p2 = NULL;

        // iterator walks through number of neighbored particles = ret_matches
        for (const auto& it : ret_matches) {

            if (!util::nearZero(it.second)) {
                p2 = &(_particles.at(it.first));

                if ((p2->getRad() + p->getRad()) > (sqrt(it.second))) {

                    // overlap
                    double overlap = (p2->getRad() + p->getRad()) - sqrt(it.second);

                    // conVec: vector from particle 1 to particle 2
                    Vector<double, 3> conVec(0., 0., 0.) ;
                    for (int i = 0; i <= 2; i++) {

                        conVec[i] = p2->getPos()[i] - p->getPos()[i];
                    }

                    Vector<double, 3> conVecNormalized(conVec) ;
                    normalize(conVecNormalized) ;

                    // Particle velocities before collision
                    Vector<double, 3> vel1bc = {p->getVel()} ;
                    Vector<double, 3> vel2bc = {p2->getVel()} ;

                    // Normal and tangential particle velocities before collision
                    Vector<double, 3> velN1(0., 0., 0.) ;
                    Vector<double, 3> velT1(0., 0., 0.) ;
                    Vector<double, 3> velN2(0., 0., 0.) ;
                    Vector<double, 3> velT2(0., 0., 0.) ;

                    // Particle velocities after collision
                    Vector<double, 3> vel1ac(0., 0., 0.) ;
                    Vector<double, 3> vel2ac(0., 0., 0.) ;

                    // Restitution coeffizient
                    double Cr = restitutionCoeff;

                    velN1 = (conVecNormalized * vel1bc) * conVecNormalized;
                    velN2 = (conVecNormalized * vel2bc) * conVecNormalized;
                    velT1 = vel1bc - velN1;
                    velT2 = vel2bc - velN2;
                    vel1ac = (Cr * p2->getMass() * ( velN2 - velN1) + (p2->getMass() * velN2) + (p->getMass() * velN1)) * (1. / (p->getMass() + p2->getMass())) ;
                    vel2ac = (Cr * p->getMass() * ( velN1 - velN2) + (p2->getMass() * velN2) + (p->getMass() * velN1)) * (1. / (p->getMass() + p2->getMass())) ;

                    double dpos[3] = {double(0), double(0), double(0) } ;

                    // Both particles are not deposited (sActivity = 3)
                    if ((p->getSActivity() != 3) && (p2->getSActivity() != 3)) {

                        for (int i = 0; i <= 2; i++) {

                            dpos[i] = conVecNormalized[i] * 0.5 * overlap ;
                            p->getPos()[i] -= dpos[i] * 1. ;
                            p->getVel()[i] = vel1ac[i] + velT1[i] ;
                            p2->getPos()[i] += dpos[i] * 1. ;
                            p2->getVel()[i] = vel2ac[i] + velT2[i] ;
                        }
                        if ((p->getSActivity() == 2) || (p->getSActivity() == 2)) {
                            p->setSActivity(2);
                            p2->setSActivity(2);
                        }
                        continue;
                    }

                    // Particle 1 is deposited (sActivity = 3) and Particle 2 is not
                    if ((p->getSActivity() != 3) && (p2->getSActivity() == 3)) {

                        for (int i = 0; i <= 2; i++) {

                            dpos[i] = conVecNormalized[i] * 0.5 * overlap ;
                            p->getPos()[i] -= 2. * dpos[i];
                            p->getVel()[i] = -1. * p->getVel()[i];
                        }
                        p->setSActivity(2) ;
                        continue;
                    }

                    // Particle 2 is deposited (sActivity = 3) and Particle 1 is not
                    if ((p->getSActivity() == 3) && (p2->getSActivity() != 3)) {

                        for (int i = 0; i <= 2; i++) {

                            dpos[i] = conVecNormalized[i] * 0.5 * overlap ;
                            p2->getPos()[i] += 2. * dpos[i];
                            p2->getVel()[i] = -1. * p2->getVel()[i];
                        }
                        p2->setSActivity(2) ;
                        continue;
                    }

                    // Both particles are deposited (sActivity = 3)
                    if ((p->getSActivity() == 3) && (p2->getSActivity() == 3)) {

                        for (int i = 0; i <= 2; i++) {

                            dpos[i] = conVecNormalized[i] * 0.5 * overlap ;
                            p->getPos()[i] -= dpos[i];
                            p2->getPos()[i] += dpos[i];
                        }
                        continue;
                    }
                }
            }
        }
    }
}

template<>
void ParticleSystem3D<double, MagneticParticle3D>::findAgglomerates()
{
    typename std::deque<MagneticParticle3D<double> >::iterator p;
    int pInt = 0;

    for (p = _particles.begin(); p != _particles.end(); ++p, ++pInt) {

        auto* p1 = &(*p);
        std::vector<std::pair<size_t, double>> ret_matches;

        getContactDetection()->getMatches(pInt, ret_matches);

        MagneticParticle3D<double>* p2 = NULL;

        for (const auto& it : ret_matches) {

            if (!util::nearZero(it.second)) {

                p2 = &(_particles.at(it.first));

                if ((p2->getRad() + p1->getRad()) > (sqrt(it.second))) {

                    // Both particles are non agglomerated
                    // A new agglomerate is formed
                    if ((p1->getAggloItr() == _Agglomerates.begin()) && (p2->getAggloItr() == _Agglomerates.begin())) {

                        std::list<MagneticParticle3D<double>*> aggloList{p1, p2} ;

                        typename std::list<MagneticParticle3D<double>*>::iterator x1;
                        typename std::list<MagneticParticle3D<double>*>::iterator x2;

                        _Agglomerates.push_back(aggloList);
                        p1->setAggloItr(_Agglomerates.end() - 1);
                        p2->setAggloItr(_Agglomerates.end() - 1);

                        for (auto x = _Agglomerates[0].begin(); x != _Agglomerates[0].end(); ++x) {

                            if (*x == p1) {
                                x1 = x;
                            }
                            if (*x == p2) {
                                x2 = x;
                            }

                        }
                        _Agglomerates[0].erase(x1);
                        _Agglomerates[0].erase(x2);

                        continue;
                    }

                    // Particle 2 is already part of an agglomerate and Particle 1 is non agglomerated
                    // Particle 1 is added to the agglomerate
                    if ((p1->getAggloItr() == _Agglomerates.begin()) && (p2->getAggloItr() != _Agglomerates.begin())) {

                        (p2->getAggloItr())->push_back(p1) ;
                        p1->setAggloItr(p2->getAggloItr()) ;
                        typename std::list<MagneticParticle3D<double>*>::iterator x1;

                        for (auto x = _Agglomerates[0].begin(); x != _Agglomerates[0].end(); ++x) {
                            if (*x == p1) {
                                x1 = x;
                            }
                        }
                        _Agglomerates[0].erase(x1);

                        continue;
                    }

                    // Particle 1 is already part of an agglomerate and Particle 2 is non agglomerated
                    // Particle 2 is added to the agglomerate
                    if ((p1->getAggloItr() != _Agglomerates.begin()) && (p2->getAggloItr() == _Agglomerates.begin())) {

                        (p1->getAggloItr())->push_back(p2) ;
                        p2->setAggloItr(p1->getAggloItr()) ;
                        typename std::list<MagneticParticle3D<double>*>::iterator x2;

                        for (auto x = _Agglomerates[0].begin(); x != _Agglomerates[0].end(); ++x) {
                            if (*x  == p2) {
                                x2 = x;
                            }
                        }

                        _Agglomerates[0].erase(x2);

                        continue;
                    }

                    // Both particles are part of diffrent agglomerates
                    // The two Agglomerates are united
                    if (((p1->getAggloItr() != _Agglomerates.begin()) && (p2->getAggloItr() != _Agglomerates.begin())) && ( p1->getAggloItr() != p2->getAggloItr())) {

                        typename std::deque<std::list<MagneticParticle3D<double>*>>::iterator x1;
                        typename std::deque<std::list<MagneticParticle3D<double>*>>::iterator x2;

                        if (p1->getAggloItr() <= p2->getAggloItr()) {
                            x2 = p2->getAggloItr() ;
                            x1 = p1->getAggloItr() ;

                        }
                        else {
                            x2 = p1->getAggloItr() ;
                            x1 = p2->getAggloItr() ;
                        }

                        x1->splice(x1->end(), *x2) ;

                        _Agglomerates.erase(x2) ;

                        for (auto anew = _Agglomerates.begin(); anew != _Agglomerates.end(); ++anew) {

                            for (auto pnew = anew->begin(); pnew != anew->end(); ++pnew) {

                                (*pnew)->setAggloItr(anew);
                            }
                        }

                        continue;

                    }
                }
            }
        }
    }
}

template<>
void ParticleSystem3D<double, MagneticParticle3D>::initAggloParticles()
{
    typename std::deque<MagneticParticle3D<double> >::iterator p;
    static int pInt = 0;

    for (p = _particles.begin() + pInt; p != _particles.end(); ++p, ++pInt) {
        p->setAggloItr(_Agglomerates.begin());
        MagneticParticle3D<double>* pPointer = &(*p);
        _Agglomerates.begin()->push_back(pPointer);
    }
}

}  // namespace olb
