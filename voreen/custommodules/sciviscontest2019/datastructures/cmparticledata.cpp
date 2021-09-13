/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2015 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "cmparticledata.h"
#include <assert.h>
#include "tgt/logmanager.h"

namespace voreen {

	int* remappingTable = nullptr;

//ParticleType getParticleType(uint16_t mask) {
//	ParticleType pType;
//	if (mask & 2) {
//		if (mask & 32) {
//			pType = BARYON_STAR;
//		}
//		if (mask & 64) {
//			pType = BARYON_WIND;
//		}
//		if (mask & 128) {
//			pType = BARYON_GAS;
//		}
//		else {
//			pType = BARYON;
//		}
//	}
//	else {
//		if (mask & 256) {
//			pType = DARKMATTER_AGN;
//		}
//		else {
//			pType = DARKMATTER;
//		}
//	}
//	return pType;
//}
/******************************************************************************
 *
 *                       CMParticleDataTimeSlice
 *
 *****************************************************************************/

    CMParticleDataTimeSlice::CMParticleDataTimeSlice(tgt::Bounds universeBounds)
            : universeBounds_(universeBounds)
    {
    }

    tgt::vec2 CMParticleDataTimeSlice::getMinMaxVelocityMagnitude(){
        int numberOfParticles       = getNumberOfParticles();
        const CMParticle* particles = startUsingParticles();
        tgt::vec2 minmax(10000000, -10000000);
        for(int i = 0; i != numberOfParticles; i++){
            float val = tgt::length(particles[i].vel);
            minmax = tgt::vec2(std::min(val, minmax.x), std::max(val, minmax.y));
        }
        finishedUsingParticles();
        return minmax;
    }

    tgt::vec2 CMParticleDataTimeSlice::getMinMaxPhi(){
        int numberOfParticles       = getNumberOfParticles();
        const CMParticle* particles = startUsingParticles();
        tgt::vec2 minmax(10000000, -10000000);
        for(int i = 0; i != numberOfParticles; i++){
            float val = particles[i].phi;
            minmax = tgt::vec2(std::min(val, minmax.x), std::max(val, minmax.y));
        }
        finishedUsingParticles();
        return minmax;
    }

	tgt::vec2 CMParticleDataTimeSlice::getMinMaxMass(){
		int numberOfParticles       = getNumberOfParticles();
		const CMParticle* particles = startUsingParticles();
		tgt::vec2 minmax(10000000, -10000000);
		for(int i = 0; i != numberOfParticles; i++){
			float val = particles[i].mass;
			minmax = tgt::vec2(std::min(val, minmax.x), std::max(val, minmax.y));
		}
		finishedUsingParticles();
		return minmax;
	}

	tgt::vec2 CMParticleDataTimeSlice::getMinMaxUU(){
		int numberOfParticles       = getNumberOfParticles();
		const CMParticle* particles = startUsingParticles();
		tgt::vec2 minmax(10000000, -10000000);
		for(int i = 0; i != numberOfParticles; i++){
			float val = particles[i].uu;
			minmax = tgt::vec2(std::min(val, minmax.x), std::max(val, minmax.y));
		}
		finishedUsingParticles();
		return minmax;
	}

	tgt::vec2 CMParticleDataTimeSlice::getMinMaxMU(){
		int numberOfParticles       = getNumberOfParticles();
		const CMParticle* particles = startUsingParticles();
		tgt::vec2 minmax(10000000, -10000000);
		for(int i = 0; i != numberOfParticles; i++){
			float val = particles[i].mu;
			minmax = tgt::vec2(std::min(val, minmax.x), std::max(val, minmax.y));
		}
		finishedUsingParticles();
		return minmax;
	}

	tgt::vec2 CMParticleDataTimeSlice::getMinMaxRho(){
		int numberOfParticles       = getNumberOfParticles();
		const CMParticle* particles = startUsingParticles();
		tgt::vec2 minmax(10000000, -10000000);
		for(int i = 0; i != numberOfParticles; i++){
			float val = particles[i].rho;
			minmax = tgt::vec2(std::min(val, minmax.x), std::max(val, minmax.y));
		}
		finishedUsingParticles();
		return minmax;
	}

	tgt::vec2 CMParticleDataTimeSlice::getMinMaxHH(){
		int numberOfParticles       = getNumberOfParticles();
		const CMParticle* particles = startUsingParticles();
		tgt::vec2 minmax(10000000, -10000000);
		for(int i = 0; i != numberOfParticles; i++){
			float val = particles[i].hh;
			minmax = tgt::vec2(std::min(val, minmax.x), std::max(val, minmax.y));
		}
		finishedUsingParticles();
		return minmax;
	}

	tgt::Bounds CMParticleDataTimeSlice::getBounds(){
		int numberOfParticles       = getNumberOfParticles();
		const CMParticle* particles = startUsingParticles();
		tgt::Bounds b;
		for(int i = 0; i != numberOfParticles; i++){
			tgt::vec3 val = particles[i].pos;
			b.addPoint(val);
		}
		finishedUsingParticles();
		return b;
	}

    tgt::Bounds CMParticleDataTimeSlice::getUniverseBounds() {
        return universeBounds_;
    }

    tgt::Bounds CMParticleDataTimeSlice::getVelocityBounds(){
        int numberOfParticles       = getNumberOfParticles();
        const CMParticle* particles = startUsingParticles();
        tgt::Bounds b;
        for(int i = 0; i != numberOfParticles; i++){
            tgt::vec3 val = particles[i].vel;
            b.addPoint(val);
        }
        finishedUsingParticles();
        return b;
    }

    tgt::mat4 CMParticleDataTimeSlice::getNormalisationTransformation() {
        //Assuming getBounds() returns bounds of the whole dataset i.e. "the universe"
        tgt::vec3 llf = getUniverseBounds().getLLF();
        //Assuming positions are currently in non comoving kpc * h0
        /*return tgt::mat4::createScale(tgt::vec3(geth0()/(1000*getTimeStep())))*
            tgt::mat4::createTranslation(-llf);*/
        return tgt::mat4::createTranslation(-llf);
    }

    int* CMParticleDataTimeSlice::buildRemappingTable(){
        int numOfParticles = getNumberOfParticles();
        //const CMParticle* particles = startUsingParticles();

        //int* remapping = new int[numOfParticles];

		if (remappingTable != nullptr) {

			remappingTable = new int[numOfParticles];
			for(int i = 0; i < numOfParticles; i++){
				/*CMParticle p = particles[i];
				assert(p.ident >= 0);
				assert(p.ident < numOfParticles);
				remapping[p.ident] = i;*/

				remappingTable[i] = i;
			}
		}

		//finishedUsingParticles();

        return remappingTable;
    }

/******************************************************************************
 *
 *                       CMParticleDataTimeSliceVector
 *
 *****************************************************************************/
    CMParticleDataTimeSliceVector::CMParticleDataTimeSliceVector(tgt::Bounds universeBounds)
            :CMParticleDataTimeSlice(universeBounds)
    {
    }
    const CMParticle* CMParticleDataTimeSliceVector::startUsingParticles(){
        return particles_.data();
    }
    int CMParticleDataTimeSliceVector::getNumberOfParticles(){
        return (int)particles_.size();
    }
    float CMParticleDataTimeSliceVector::getTimeStep(){
        return timeStep_;
    }

    const int* CMParticleDataTimeSliceVector::getRemappingTable(){
        int numOfParticles = getNumberOfParticles();
        if (remappingTable_.size() < numOfParticles){
            remappingTable_.resize(numOfParticles);
            for(int i = 0; i != numOfParticles; i++){
                remappingTable_[i] = i;
            }
        }
        return remappingTable_.data();
    }

    float CMParticleDataTimeSliceVector::geth0() {
        //std::cout << "Fix geth0 in particleslicevector!" << std::endl;
        //return 0.688062;
        return h0_;
    }

/******************************************************************************
 *
 *                                CMParticleData
 *
 *****************************************************************************/
    const std::string CMParticleData::loggerCat_("voreen.CMParticleData");

    bool timeSliceIsEarlier(CMParticleDataTimeSlice* t1, CMParticleDataTimeSlice* t2) {
        return t1->getTimeStep() < t2->getTimeStep();
    }

    CMParticleData::CMParticleData(const std::vector<CMParticleDataTimeSlice*>& timeSlices)
            : timeSlices_(timeSlices)
            , isFiltered_(false)
    {
        std::sort(timeSlices_.begin(), timeSlices_.end(), timeSliceIsEarlier);
        int numOfParticles = timeSlices[0]->getNumberOfParticles();
        /*for(auto slice: timeSlices){
            tgtAssert(slice->getNumberOfParticles() == numOfParticles, "CMParticleData: All timeslices need same number of particles");
        }*/
        LDEBUG("alloc pd");
    }
    CMParticleData::~CMParticleData() {
        for(auto it : timeSlices_)
            delete it;
        timeSlices_.clear();
        LDEBUG("free pd");
    }
    const std::vector<CMParticleDataTimeSlice*>& CMParticleData::particleDataTimeSlices() const {
        return timeSlices_;
    }

    CMParticleDataTimeSlice* CMParticleData::sliceAtTimeStep(float a) const {
        if(timeSlices_.empty()) {
            return nullptr;
        }
        CMParticleDataTimeSlice* lastData = timeSlices_[(int)roundf(a)];
        //We know that timeSlices_ is sorted, so we can just go from left to right:
        /*for(size_t i=1; i < timeSlices_.size(); ++i) {
            if(lastData->getTimeStep() + timeSlices_[i]->getTimeStep() > 2*a) {
                if(lastData == nullptr) {
                    return timeSlices_[i];
                } else {
                    return lastData;
                }
            } else {
                lastData = timeSlices_[i];
            }
        }*/
        return lastData;
    }


    const char* CMParticleData::getEnabledState() const{
        // FIXME!
        CMParticleData * pd = const_cast<CMParticleData*>(this);
        bool filtered = isFiltered_;
        const char* c = pd->getEnabledStateMutable();
        pd->isFiltered_ = filtered;
        return c;
    }

    char* CMParticleData::getEnabledStateMutable(){
        int numOfParticles = timeSlices_[0]->getNumberOfParticles();
        if (enabledState_.size() != numOfParticles){
            enabledState_.resize(numOfParticles);
            for(int i = 0; i != numOfParticles; i++){
                enabledState_[i] = 1;
            }
        }
        isFiltered_ = true;
        return enabledState_.data();
    }

    namespace{

// THIS IS A HORRIBLE IDEA
        class PasstroughTimeSlice : public CMParticleDataTimeSlice{
        public:
            PasstroughTimeSlice(CMParticleDataTimeSlice* timeSlice);

            virtual const CMParticle* startUsingParticles() override;
            virtual void              finishedUsingParticles() override;
            virtual int               getNumberOfParticles() override;
            virtual float             getTimeStep() override;
            virtual const int*        getRemappingTable() override;
            virtual float             geth0() override;
            virtual tgt::mat4         getNormalisationTransformation() override;
        private:
            CMParticleDataTimeSlice* timeSlice_;
        };

        tgt::mat4 PasstroughTimeSlice::getNormalisationTransformation(){
            return timeSlice_->getNormalisationTransformation();
        }

        PasstroughTimeSlice::PasstroughTimeSlice(CMParticleDataTimeSlice* timeSlice)
                : CMParticleDataTimeSlice(timeSlice->getUniverseBounds())
                , timeSlice_(timeSlice)
        {
        }


        const CMParticle* PasstroughTimeSlice::startUsingParticles(){
            return timeSlice_->startUsingParticles();
        }

        void PasstroughTimeSlice::finishedUsingParticles(){
            timeSlice_->finishedUsingParticles();
        }

        int PasstroughTimeSlice::getNumberOfParticles(){
            return timeSlice_->getNumberOfParticles();
        }

        float PasstroughTimeSlice::getTimeStep(){
            return timeSlice_->getTimeStep();
        }

        const int* PasstroughTimeSlice::getRemappingTable(){
            return timeSlice_->getRemappingTable();;
        }
        float PasstroughTimeSlice::geth0() {
            return timeSlice_->geth0();
        }

    }


    bool CMParticleData::isFiltered() const{
        return isFiltered_;
    }

    CMParticleData* CMParticleData::cloneWithParticles() const{
        std::vector<CMParticleDataTimeSlice*> timeSlices;
        for(auto slice:timeSlices_){
            timeSlices.push_back(new PasstroughTimeSlice(slice));
        }
        CMParticleData* particleData = new CMParticleData(timeSlices);

        int numOfParticles = timeSlices[0]->getNumberOfParticles();

        if (isFiltered()){
            const char* oldEnabled = getEnabledState();
            char* newEnabled = particleData->getEnabledStateMutable();
            for(int i = 0; i != numOfParticles; i++){
                newEnabled[i] = oldEnabled[i];
            }
        }
        return particleData;
    }
} // namespace
