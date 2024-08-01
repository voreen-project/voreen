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

// header
#include "cosmologyparticlefilter.h"

// addtl. headers
#include "tgt/filesystem.h"

namespace voreen {

    using tgt::vec3;
    using tgt::MouseEvent;
    using std::vector;


    //inport_.onNewData(MemberFunctionCallback<CosmologyParticleFilter>(this, &CosmologyParticleFilter::process));

    CosmologyParticleFilter::CosmologyParticleFilter()
        : Processor()
        , outport_(Port::OUTPORT, "particlehandle.output", "Particle Data Output")
        , inport_(Port::INPORT, "particlehandle.input", "Particle Data Input")
        , timeStep_("timeStep", "Filter at Time Step", 0.0f, 0.0f, 624.0f)
        , linkedTimeStep_("linkedTimeStep", "Linked Time Step", 0.0f, 0.0f, 624.0f)
        , useLinkedTime_("uselinked", "Use Linked Time Step", false)
        , currentfilter_("currentfilter", "Current Filter", 0, 0, 9)
        , enabled_("enabled", "Enable", false)
        , anyTimeFilter_("anyTime", "Any time filter", false)
        , invert_("invert", "Invert", false)
        , filtertype_("filtertype", "Filter Type")
        , filterconnection_("filterconnection", "Connect filters by")
        , intervalfilter_("intervalfilter", "Interval type")
        , pos_("position", "Position", tgt::vec3(0.0f))
        , normal_("normal", "Normal vector", tgt::vec3(0.0f))
        , radius_("radius", "Radius")
        , bounds_("bounds", "Bounds of Bounding Box")
        , interval_("interval", "Interval", tgt::vec2(0.0f))
        , evalButton_("evalButton", "Evaluate filter")
        , doFilterOnNextEvaluation_(true)
		//, typeSelection_("typeselection", "Select Types", std::vector<Property*>())
    {
        addPort(outport_);
        addPort(inport_);

        addProperty(timeStep_);
        addProperty(linkedTimeStep_);
		linkedTimeStep_.setVisibleFlag(false);
        addProperty(useLinkedTime_);
		ON_CHANGE_LAMBDA(useLinkedTime_, [this] {timeStep_.setVisibleFlag(!useLinkedTime_.get()); });
		ON_CHANGE_LAMBDA(linkedTimeStep_, [this] {if (!doFilterOnNextEvaluation_) { doFilterOnNextEvaluation_ = useLinkedTime_.get(); } });


        addProperty(currentfilter_);
        addProperty(filterconnection_);
        addProperty(enabled_);
        addProperty(anyTimeFilter_);
        addProperty(invert_);
        addProperty(filtertype_);
        addProperty(pos_);
        addProperty(normal_);
        addProperty(radius_);
        addProperty(bounds_);
        addProperty(intervalfilter_);
        addProperty(interval_);
		//addProperty(typeSelection_);


        filtertype_.addOption("trivial", "Trivial", trivial);
		filtertype_.addOption("type", "Particle Type", particletype);
        filtertype_.addOption("plane", "Plane", plane);
        filtertype_.addOption("interval", "Interval", interval);
        filtertype_.addOption("box", "Box", box);
        filtertype_.addOption("sphere", "Sphere", sphere);

		filterconnection_.addOption("and", "AND", fc_and);
		filterconnection_.addOption("or", "OR", fc_or);

        intervalfilter_.addOption("posx", "X Position", posx);
        intervalfilter_.addOption("posy", "Y Position", posy);
        intervalfilter_.addOption("posz", "Z Position", posz);
        intervalfilter_.addOption("velx", "X Velocity", velx);
        intervalfilter_.addOption("vely", "Y Velocity", vely);
        intervalfilter_.addOption("velz", "Z Velocity", velz);
        intervalfilter_.addOption("tvel", "Total Velocity", tvel);
        intervalfilter_.addOption("phi", "Phi value", phi);
        intervalfilter_.addOption("mass", "Mass", mass);
        intervalfilter_.addOption("uu", "Internal energy", uu);
        intervalfilter_.addOption("mu", "Molecular weight", mu);
        intervalfilter_.addOption("rho", "Density", rho);
        intervalfilter_.addOption("hh", "Smoothing length", hh);

		allDarkMatter_ = BoolProperty("darkmatter", "All dark matter", false);
		allBaryons_ = BoolProperty("baryons", "All baryons", false);
		starBaryons_ = BoolProperty("starbaryons", "Star Particles", false);
		windBaryons_ = BoolProperty("windbaryons", "Wind Particles", false);
		gasBaryons_ = BoolProperty("gasbaryons", "Gas Particles", false);
		agn_ = BoolProperty("agn", "Active galactic nuclei", false);

		allBaryons_.setGroupName("Select Types:");
		allBaryons_.setGroupID("type");
		allDarkMatter_.setGroupID("type");
		starBaryons_.setGroupID("type");
		windBaryons_.setGroupID("type");
		gasBaryons_.setGroupID("type");
		agn_.setGroupID("type");
		/*typeSelection_.addProperty(allDarkMatter_);
		typeSelection_.addProperty(allBaryons_);
		typeSelection_.addProperty(starBaryons_);
		typeSelection_.addProperty(windBaryons_);
		typeSelection_.addProperty(gasBaryons_);
		typeSelection_.addProperty(agn_);*/

		addProperty(allBaryons_);
		addProperty(allDarkMatter_);
		addProperty(starBaryons_);
		addProperty(windBaryons_);
		addProperty(gasBaryons_);
		addProperty(agn_);

        ON_CHANGE(currentfilter_, CosmologyParticleFilter, switchCurrentFilter);
        ON_CHANGE(currentfilter_, CosmologyParticleFilter, switchFilterType);
        ON_CHANGE(filtertype_, CosmologyParticleFilter, switchFilterType)
        ON_CHANGE(timeStep_, CosmologyParticleFilter, onEvalButtonPressed);
        //ON_CHANGE(autoEval_, CosmologyParticleFilter, onAutoEvalChanged);
        ON_CHANGE(evalButton_, CosmologyParticleFilter, onEvalButtonPressed);

		blockTypeCallBacks(true);
		ON_CHANGE_LAMBDA(allDarkMatter_, [this] {onTypeChanged(alldarkmatter); });
		ON_CHANGE_LAMBDA(allBaryons_, [this] {onTypeChanged(allbaryons); });
		ON_CHANGE_LAMBDA(starBaryons_, [this] {onTypeChanged(starparticle); });
		ON_CHANGE_LAMBDA(windBaryons_, [this] {onTypeChanged(windparticle); });
		ON_CHANGE_LAMBDA(gasBaryons_, [this] {onTypeChanged(gasparticle); });
		ON_CHANGE_LAMBDA(agn_, [this] {onTypeChanged(agn); });
		blockTypeCallBacks(false);
		//agn_.onChange([this] {this.onTypeChanged(1); })

        //inport_.onNewData(MemberFunctionCallback<CosmologyParticleFilter>(this, &CosmologyParticleFilter::onEvalButtonPressed));
        
		pos_.setVisibleFlag(false);
		normal_.setVisibleFlag(false);
		radius_.setVisibleFlag(false);
		bounds_.setVisibleFlag(false);
		intervalfilter_.setVisibleFlag(false);
		interval_.setVisibleFlag(false);
		setPropertyGroupVisible("type", false);

        addProperty(evalButton_);
    }
    CosmologyParticleFilter::~CosmologyParticleFilter(){
    }
	
    void CosmologyParticleFilter::initialize() {
        // call superclass function first
        Processor::initialize();
        // Class is initialized with a single, empty filter
        CosmologyParticleFilter::addFilter();
        currentfilter_.set(0);
        CosmologyParticleFilter::popFilter();
    }

    Processor* CosmologyParticleFilter::create() const{
        return new CosmologyParticleFilter;
    }

    void CosmologyParticleFilter::deinitialize() {
        // call superclass function last
        filters_.clear();
        Processor::deinitialize();
    }

    bool CosmologyParticleFilter::isReady() const {
        return outport_.isReady();
    }
    
    void CosmologyParticleFilter::serialize(Serializer& s) const {
    	// call superclass function first
        Processor::serialize(s);
        s.serialize("filters", filters_);
    }


    void CosmologyParticleFilter::deserialize(Deserializer& s){
    	// call superclass function first
		blockTypeCallBacks(true);

        Processor::deserialize(s);
		blockTypeCallBacks(false);

        s.optionalDeserialize("filters", filters_, std::vector<PartFilter>());

		switchFilterType();
    }


    void CosmologyParticleFilter::addFilter(){
        PartFilter p = PartFilter();
        filters_.push_back(p);
    }

    void CosmologyParticleFilter::removeLastFilter(){
        filters_.pop_back();

    }

    void CosmologyParticleFilter::popFilter(){
        int i = currentfilter_.get();

        invert_.set(filters_[i].invert_);
        enabled_.set(filters_[i].enabled_);
		anyTimeFilter_.set(filters_[i].anyTime_);
        filtertype_.selectByValue(filters_[i].filtertype_);
        intervalfilter_.selectByValue(filters_[i].intervalfilter_);
        pos_.set(filters_[i].pos_);
        normal_.set(filters_[i].normal_);
        radius_.set(filters_[i].radius_);
        bounds_.set(filters_[i].bounds_);
        interval_.set(filters_[i].interval_);

		blockTypeCallBacks(true);
		allDarkMatter_.set(filters_[i].allDarkMatter_);
		allBaryons_.set(filters_[i].allBaryons_);
		starBaryons_.set(filters_[i].starBaryons_);
		windBaryons_.set(filters_[i].windBaryons_);
		gasBaryons_.set(filters_[i].gasBaryons_);
		agn_.set(filters_[i].agn_);
		blockTypeCallBacks(false);

    }

    void CosmologyParticleFilter::pushFilter(){
        int i = currentfilter_.get();

        filters_[i].invert_ = invert_.get();
        filters_[i].enabled_ = enabled_.get();
		filters_[i].anyTime_ = anyTimeFilter_.get();
        filters_[i].filtertype_ = filtertype_.getValue();
        filters_[i].intervalfilter_ = intervalfilter_.getValue();
        filters_[i].pos_ = pos_.get();
        filters_[i].normal_ = normal_.get();
        filters_[i].radius_ = radius_.get();
        filters_[i].bounds_ = bounds_.get();
        filters_[i].interval_ = interval_.get();
		filters_[i].allDarkMatter_ = allDarkMatter_.get();
		filters_[i].allBaryons_ = allBaryons_.get();
		filters_[i].starBaryons_ = starBaryons_.get();
		filters_[i].windBaryons_ = windBaryons_.get();
		filters_[i].gasBaryons_ = gasBaryons_.get();
		filters_[i].agn_ = agn_.get();

    }
    void CosmologyParticleFilter::switchCurrentFilter(){
        int i = currentfilter_.get();
        if (i < 0){
            currentfilter_.set(0);
            return;
        }
            while (filters_.size() <= i){
                CosmologyParticleFilter::addFilter();
            }
            popFilter();

    }



    void CosmologyParticleFilter::switchFilterType(){
        if (filters_.size()){
            switch (filtertype_.getValue()){
            case (trivial) :
                pos_.setVisibleFlag(false);
                normal_.setVisibleFlag(false);
                radius_.setVisibleFlag(false);
                bounds_.setVisibleFlag(false);
                intervalfilter_.setVisibleFlag(false);
                interval_.setVisibleFlag(false);
				//typeSelection_.setVisibleFlag(false);
				setPropertyGroupVisible("type", false);
                break;

            case (plane) :
                pos_.setVisibleFlag(true);
                normal_.setVisibleFlag(true);
                radius_.setVisibleFlag(false);
                bounds_.setVisibleFlag(false);
                intervalfilter_.setVisibleFlag(false);
                interval_.setVisibleFlag(false);
				//typeSelection_.setVisibleFlag(false);
				setPropertyGroupVisible("type", false);
                break;

            case (interval) :
                pos_.setVisibleFlag(false);
                normal_.setVisibleFlag(false);
                radius_.setVisibleFlag(false);
                bounds_.setVisibleFlag(false);
                intervalfilter_.setVisibleFlag(true);
                interval_.setVisibleFlag(true);
				//typeSelection_.setVisibleFlag(false);
				setPropertyGroupVisible("type", false);
                break;

            case (box) :
                pos_.setVisibleFlag(false);
                normal_.setVisibleFlag(false);
                radius_.setVisibleFlag(false);
                bounds_.setVisibleFlag(true);
                intervalfilter_.setVisibleFlag(false);
                interval_.setVisibleFlag(false);
				//typeSelection_.setVisibleFlag(false);
				setPropertyGroupVisible("type", false);
                break;

            case (sphere) :
                pos_.setVisibleFlag(true);
                normal_.setVisibleFlag(false);
                radius_.setVisibleFlag(true);
                bounds_.setVisibleFlag(false);
                intervalfilter_.setVisibleFlag(false);
                interval_.setVisibleFlag(false);
				//typeSelection_.setVisibleFlag(false);
				setPropertyGroupVisible("type", false);
                break;

			case (particletype):
				pos_.setVisibleFlag(false);
				normal_.setVisibleFlag(false);
				radius_.setVisibleFlag(false);
				bounds_.setVisibleFlag(false);
				intervalfilter_.setVisibleFlag(false);
				interval_.setVisibleFlag(false);
				//typeSelection_.setVisibleFlag(true);
				setPropertyGroupVisible("type", true);
				break;
            }
        }

    }


    void CosmologyParticleFilter::process() {
        if (!doFilterOnNextEvaluation_)
            return;
        const CMParticleData* oldParticleData = inport_.getData();
        if (!oldParticleData)
            return;

		// Access the particle data and create a copy
        CMParticleData*           newParticleData = oldParticleData->cloneWithParticles();

		// Access the current time step

        

        // Access the Array containing the information which particles are active,
        // as Booleans reaching from 0 to numOfParticles
        char*                     enabledState = newParticleData->getEnabledStateMutable();

		int numOfParticles = newParticleData->sliceAtTimeStep(0.0f)->getNumberOfParticles();

		bool reset = true;

		if (filterconnection_.getValue() != fc_and ) {
			memset(enabledState, 0, numOfParticles);
		}

        // Cycle through filters
        for (int k = 0; k < filters_.size(); k++){
            if (filters_[k].enabled_) {
				if (filters_[k].anyTime_) {

					char* tmpEnabled = new char[numOfParticles];
					memset(tmpEnabled, 0, numOfParticles);

					for (int step = 0; step < 625; step++) {
						CMParticleDataTimeSlice* newTimeSlice = newParticleData->sliceAtTimeStep((float)step);
						doFilter(newTimeSlice, tmpEnabled, k, false);
					}

					for (int i = 0; i < numOfParticles; i++) {
						if (filterconnection_.getValue() == fc_and) {
							enabledState[i] = enabledState[i] && tmpEnabled[i];
						}
						else {
							enabledState[i] = enabledState[i] || tmpEnabled[i];
						}
					}
					//memcpy(enabledState, tmpEnabled, numOfParticles);
					delete[] tmpEnabled;
 				}
				else {
					/*if (filterconnection_.getValue() != fc_and &&reset) {
						memset(enabledState, 0, numOfParticles);
					}*/
					float timeStep = (useLinkedTime_.get()) ? linkedTimeStep_.get() : timeStep_.get();
					CMParticleDataTimeSlice * newTimeSlice = newParticleData->sliceAtTimeStep(timeStep);
					doFilter(newTimeSlice, enabledState, k, filterconnection_.getValue() == fc_and);

					reset = false;
				}
                
            }
        }
        

        outport_.setData(newParticleData);

        doFilterOnNextEvaluation_ = false;
    }

	void CosmologyParticleFilter::doFilter(CMParticleDataTimeSlice* particleData, char* enabledState, int k, bool connectByAnd){
		// Prepare to convert relative values into coordinates, based on the current timestep
		tgt::mat4 M;
		M = CMinvert(particleData->getNormalisationTransformation());
		// Access current absolute values
		tgt::Bounds universeBounds_ = particleData->getUniverseBounds();
		tgt::Bounds vel_ = tgt::Bounds();
		tgt::vec2 mmvel_ = tgt::vec2();
		tgt::vec2 phi_ = tgt::vec2();
		tgt::vec2 mass_ = tgt::vec2();
		tgt::vec2 uu_ = tgt::vec2();
		tgt::vec2 mu_ = tgt::vec2();
		tgt::vec2 rho_ = tgt::vec2();
		tgt::vec2 hh_ = tgt::vec2();

		if (filters_[k].filtertype_ == interval) {
			switch (filters_[k].intervalfilter_) {
				// compute true values of interval bounds

			case (velx):
			case (vely):
			case (velz):
				vel_ = particleData->getVelocityBounds();
				break;
			case (tvel):
				mmvel_ = particleData->getMinMaxVelocityMagnitude();
				break;
			case (phi):
				phi_ = particleData->getMinMaxPhi();
				break;
			case (mass):
				mass_ = particleData->getMinMaxMass();
				break;
			case (uu):
				uu_ = particleData->getMinMaxUU();
				break;
			case (mu):
				mu_ = particleData->getMinMaxMU();
				break;
			case (rho):
				rho_ = particleData->getMinMaxRho();
				break;
			case (hh):
				hh_ = particleData->getMinMaxHH();
				break;
			}
		}
		const CMParticle* particles = particleData->startUsingParticles();

		int numOfParticles = particleData->getNumberOfParticles();

		// Filter particles
#pragma omp parallel for
		for (int i = 0; i < numOfParticles; i++) {
			CMParticle p = particles[i];
			bool visitParticle = (connectByAnd) ? enabledState[p.ident] : true;

			if (visitParticle) {
				bool b_;		// filter evaluation variable
				switch (filters_[k].filtertype_) {
				case (trivial): b_ = false;
					break;

				case (particletype):
				{
					uint16_t mask = 0;
					if (filters_[k].allBaryons_ && filters_[k].allDarkMatter_) {
						b_ = true;
					}
					else if (filters_[k].allDarkMatter_) {
						mask = (filters_[k].starBaryons_ ? 32 : 0) + (filters_[k].windBaryons_ ? 64 : 0) + (filters_[k].gasBaryons_ ? 128 : 0);

						b_ = (p.mask & mask || !(p.mask & 2)) != 0;
					}
					else if (filters_[k].allBaryons_) {
						mask = 2;

						b_ = (p.mask & mask || (p.mask & (filters_[k].agn_ ? 256 : 0))) != 0;
					}
					else {
						mask = ((filters_[k].starBaryons_ ? 32 : 0) + (filters_[k].windBaryons_ ? 64 : 0) + (filters_[k].gasBaryons_ ? 128 : 0) + (filters_[k].agn_ ? 256 : 0));

						b_ = (p.mask & mask) != 0;
					}
				}
				break;

				case (plane):
					filters_[k].truepos_ = universeBounds_.getLLF() * (tgt::vec3(0.0f) - filters_[k].pos_) + filters_[k].pos_ * universeBounds_.getURB();
					b_ = (tgt::dot((p.pos - filters_[k].truepos_), filters_[k].normal_) <= 0);
					break;

				case (interval): {
					float f;
					switch (filters_[k].intervalfilter_) {
						// compute true values of interval bounds
					case (posx):
						f = p.pos[0];
						filters_[k].trueinterval_[0] = universeBounds_.getLLF()[0] * (1.0f - filters_[k].interval_[0]) + filters_[k].interval_[0] * universeBounds_.getURB()[0];
						filters_[k].trueinterval_[1] = universeBounds_.getLLF()[0] * (1.0f - filters_[k].interval_[1]) + filters_[k].interval_[1] * universeBounds_.getURB()[0];
						break;
					case (posy):
						f = p.pos[1];
						filters_[k].trueinterval_[0] = universeBounds_.getLLF()[1] * (1.0f - filters_[k].interval_[0]) + filters_[k].interval_[0] * universeBounds_.getURB()[1];
						filters_[k].trueinterval_[1] = universeBounds_.getLLF()[1] * (1.0f - filters_[k].interval_[1]) + filters_[k].interval_[1] * universeBounds_.getURB()[1];
						break;
					case (posz):
						f = p.pos[2];
						filters_[k].trueinterval_[0] = universeBounds_.getLLF()[2] * (1.0f - filters_[k].interval_[0]) + filters_[k].interval_[0] * universeBounds_.getURB()[2];
						filters_[k].trueinterval_[1] = universeBounds_.getLLF()[2] * (1.0f - filters_[k].interval_[1]) + filters_[k].interval_[1] * universeBounds_.getURB()[2];
						break;
					case (velx):
						f = p.vel[0];
						filters_[k].trueinterval_[0] = vel_.getLLF()[0] * (1.0f - filters_[k].interval_[0]) + filters_[k].interval_[0] * vel_.getURB()[0];
						filters_[k].trueinterval_[1] = vel_.getLLF()[0] * (1.0f - filters_[k].interval_[1]) + filters_[k].interval_[1] * vel_.getURB()[0];
						break;
					case (vely):
						f = p.vel[1];
						filters_[k].trueinterval_[0] = vel_.getLLF()[1] * (1.0f - filters_[k].interval_[0]) + filters_[k].interval_[0] * vel_.getURB()[1];
						filters_[k].trueinterval_[1] = vel_.getLLF()[1] * (1.0f - filters_[k].interval_[1]) + filters_[k].interval_[1] * vel_.getURB()[1];
						break;
					case (velz):
						f = p.vel[2];
						filters_[k].trueinterval_[0] = vel_.getLLF()[2] * (1.0f - filters_[k].interval_[0]) + filters_[k].interval_[0] * vel_.getURB()[2];
						filters_[k].trueinterval_[1] = vel_.getLLF()[2] * (1.0f - filters_[k].interval_[1]) + filters_[k].interval_[1] * vel_.getURB()[2];
						break;
					case (tvel):
						f = tgt::length(p.vel);
						filters_[k].trueinterval_[0] = mmvel_[0] * (1.0f - filters_[k].interval_[0]) + filters_[k].interval_[0] * mmvel_[1];
						filters_[k].trueinterval_[1] = mmvel_[0] * (1.0f - filters_[k].interval_[1]) + filters_[k].interval_[1] * mmvel_[1];
						break;
					case (phi):
						f = p.phi;
						filters_[k].trueinterval_[0] = phi_[0] * (1.0f - filters_[k].interval_[0]) + filters_[k].interval_[0] * phi_[1];
						filters_[k].trueinterval_[1] = phi_[0] * (1.0f - filters_[k].interval_[1]) + filters_[k].interval_[1] * phi_[1];
						break;
					case (mass):
						f = p.mass;
						filters_[k].trueinterval_[0] = mass_[0] * (1.0f - filters_[k].interval_[0]) + filters_[k].interval_[0] * mass_[1];
						filters_[k].trueinterval_[1] = mass_[0] * (1.0f - filters_[k].interval_[1]) + filters_[k].interval_[1] * mass_[1];
						break;
					case (uu):
						f = p.uu;
						filters_[k].trueinterval_[0] = uu_[0] * (1.0f - filters_[k].interval_[0]) + filters_[k].interval_[0] * uu_[1];
						filters_[k].trueinterval_[1] = uu_[0] * (1.0f - filters_[k].interval_[1]) + filters_[k].interval_[1] * uu_[1];
						break;
					case (mu):
						f = p.mu;
						filters_[k].trueinterval_[0] = mu_[0] * (1.0f - filters_[k].interval_[0]) + filters_[k].interval_[0] * mu_[1];
						filters_[k].trueinterval_[1] = mu_[0] * (1.0f - filters_[k].interval_[1]) + filters_[k].interval_[1] * mu_[1];
						break;
					case (rho):
						f = p.rho;
						filters_[k].trueinterval_[0] = rho_[0] * (1.0f - filters_[k].interval_[0]) + filters_[k].interval_[0] * rho_[1];
						filters_[k].trueinterval_[1] = rho_[0] * (1.0f - filters_[k].interval_[1]) + filters_[k].interval_[1] * rho_[1];
						break;
					case (hh):
						f = p.hh;
						filters_[k].trueinterval_[0] = hh_[0] * (1.0f - filters_[k].interval_[0]) + filters_[k].interval_[0] * hh_[1];
						filters_[k].trueinterval_[1] = hh_[0] * (1.0f - filters_[k].interval_[1]) + filters_[k].interval_[1] * hh_[1];
						break;
					}

					float upper = std::max(filters_[k].trueinterval_[0], filters_[k].trueinterval_[1]);
					float lower = std::min(filters_[k].trueinterval_[0], filters_[k].trueinterval_[1]);
					b_ = ((lower <= f) & (f <= upper));
				}
								 break;

				case (box):
					filters_[k].truebounds_ = CMtransformBounds(M, filters_[k].bounds_);
					b_ = filters_[k].truebounds_.containsPoint(p.pos);
					break;
				case (sphere):
					filters_[k].truepos_ = universeBounds_.getLLF() * (tgt::vec3(0.0f) - filters_[k].pos_) + filters_[k].pos_ * universeBounds_.getURB();
					filters_[k].trueradius_ = filters_[k].radius_ * tgt::length(particleData->getUniverseBounds().diagonal()) * 0.5f;
					b_ = (tgt::distance(p.pos, filters_[k].truepos_) <= filters_[k].trueradius_);
					break;
				}
				// apply filter evaluation result
				if (connectByAnd) {
					enabledState[p.ident] = (b_ ^ filters_[k].invert_);
				} 
				else {
					enabledState[p.ident] = enabledState[p.ident] || (b_);// ^ filters_[k].invert_);
				}
			}
		}
		particleData->finishedUsingParticles();
	}

    void CosmologyParticleFilter::onEvalButtonPressed(){
//        filters_.back().enabled_ = enabled_.get();
//        filters_.back().invert_ = invert_.get();
        if (filters_.size()){
            pushFilter();
        }
        doFilterOnNextEvaluation_ = true;
        invalidate();
    }

	void CosmologyParticleFilter::onTypeChanged(TYPECHANGED box) {
		//LWARNING("Set type . ");
		switch (box) {
		case alldarkmatter: {
				bool setTo = allDarkMatter_.get();

				if (setTo) {
					agn_.set(false);
				}
			}
			break;
		case allbaryons: {
				bool setTo = allBaryons_.get();

				if (setTo) {
					starBaryons_.set(false);
					windBaryons_.set(false);
					gasBaryons_.set(false);
				}
			}
			break;
		case starparticle: {
				bool setTo = starBaryons_.get();

				if (setTo) {
					allBaryons_.set(false);
				}
			}
			break;
		case windparticle: {
				bool setTo = windBaryons_.get();

				if (setTo) {
					allBaryons_.set(false);
				}
			}
			break;
		case gasparticle: {
				bool setTo = gasBaryons_.get();

				if (setTo) {
					allBaryons_.set(false);
				}
			}
			break;
		case agn: {
				bool setTo = agn_.get();

				if (setTo) {
					allDarkMatter_.set(false);
				}
			}
			break;
		}
	}

	//custom serialization and deserialization
    void CosmologyParticleFilter::PartFilter::serialize(Serializer& s) const{
        s.serialize("enabled", enabled_);
        s.serialize("anyTime", anyTime_);
        int filtertype = static_cast<int>(filtertype_);
        s.serialize("filtertype", filtertype);
        int intervalfilter = static_cast<int>(intervalfilter_);
        s.serialize("intervalfilter", intervalfilter);
        s.serialize("pos", pos_);
        s.serialize("truepos", truepos_);
        s.serialize("normal", normal_);
        s.serialize("radius", radius_);
        s.serialize("trueradius", trueradius_);
        s.serialize("bounds", bounds_);
        s.serialize("truebounds", truebounds_);
        s.serialize("interval", interval_);
        s.serialize("trueinterval", trueinterval_);

		s.serialize("darkmatter", allDarkMatter_);
		s.serialize("baryons", allBaryons_);
		s.serialize("stars", starBaryons_);
		s.serialize("wind", windBaryons_);
		s.serialize("gas", gasBaryons_);
		s.serialize("agn", agn_);

    }
    void CosmologyParticleFilter::PartFilter::deserialize(Deserializer& s){
        s.optionalDeserialize("enabled", enabled_, false);
        s.optionalDeserialize("anyTime", anyTime_, false);
        int filtertype;
        s.optionalDeserialize("filtertype", filtertype, 0);
        filtertype_ = static_cast<FILTERTYPE>(filtertype);
        int intervalfilter;
        s.optionalDeserialize("intervalfilter", intervalfilter, 0);
        intervalfilter_ = static_cast<INTERVALFILTER>(intervalfilter);
        s.optionalDeserialize("pos", pos_, tgt::vec3::zero);
        s.optionalDeserialize("truepos", truepos_, tgt::vec3::zero);
        s.optionalDeserialize("normal", normal_, tgt::vec3::zero);
        s.optionalDeserialize("radius", radius_, 1.0f);
        s.optionalDeserialize("trueradius", trueradius_, 1.0f);
        s.optionalDeserialize("bounds", bounds_, tgt::Bounds());
        s.optionalDeserialize("truebounds", truebounds_, tgt::Bounds());
        s.optionalDeserialize("interval", interval_, tgt::vec2(0, 1));
        s.optionalDeserialize("trueinterval", trueinterval_, tgt::vec2(0, 1));

		s.optionalDeserialize("darkmatter", allDarkMatter_, false);
		s.optionalDeserialize("baryons", allBaryons_, false);
		s.optionalDeserialize("stars", starBaryons_, false);
		s.optionalDeserialize("wind", windBaryons_, false);
		s.optionalDeserialize("gas", gasBaryons_, false);
		s.optionalDeserialize("agn", agn_, false);
    }

}   // end of namespace
