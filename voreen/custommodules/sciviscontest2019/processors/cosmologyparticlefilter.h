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

#ifndef VRN_COSMOLOGYPARTICLEFILTER_H
#define VRN_COSMOLOGYPARTICLEFILTER_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/propertyvector.h"

#include "../ports/cmparticleport.h"

#include "../datastructures/cmparticledata.h"
#include "../utils/cmmath.h"

namespace voreen{
    /**
    * A throughput processor for cosmology particle data. It takes data from a cosmology 
    * data source and manipulates the enabledState flag of the particles according to a 
	* configurable set of binary restrictions. These filters can be manipulated and
	* activated independently and are persistent for each workspace.	
    */
    class CosmologyParticleFilter : public Processor{
    public:
        CosmologyParticleFilter();
        virtual ~CosmologyParticleFilter();
        virtual Processor* create() const;

        virtual std::string getClassName() const    { return "CosmologyParticleFilter"; }
        virtual std::string getCategory() const     { return "Viscontest2019"; }
        virtual void setDescriptions()              { setDescription("Cosmology Particle Filter"); }
        virtual CodeState   getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }
        virtual bool isReady() const;               

        virtual void serialize(Serializer& s) const override;
        virtual void deserialize(Deserializer& s) override;

    protected:
        virtual void initialize() override;
        virtual void deinitialize();

        virtual void process();
        
	/**
	 *	Implementation of the PartFilter object which contains all information 
	 *	pertaining to a specific filter.
	 */
    private:
        enum FILTERTYPE { trivial, plane, interval, box, sphere, particletype};
        enum FILTERCONNECTION { fc_and, fc_or};
        enum INTERVALFILTER { posx, posy, posz, velx, vely, velz, tvel, phi, mass, uu, mu, rho, hh };
		enum TYPECHANGED { alldarkmatter, allbaryons, starparticle, windparticle, gasparticle, agn };

        struct PartFilter :public Serializable{
        public:
            PartFilter()
                : invert_(false)
                , enabled_(false)
                , anyTime_(false)
                , filtertype_(trivial)
                , intervalfilter_(phi)
                , pos_(tgt::vec3::zero)
                , normal_(tgt::vec3::zero)
                , radius_(0.0f)
                //, bounds_(tgt::Bounds)
                , interval_(tgt::vec2::zero)
				, allDarkMatter_(false)
				, allBaryons_(false)
				, starBaryons_(false)
				, windBaryons_(false)
				, gasBaryons_(false)
				, agn_(false)
            {
            }
            
     /**
      *	Initialization of the "currently active" filter, for runtime user interaction
      */

            bool invert_;
            bool enabled_;
            bool anyTime_;
            FILTERTYPE filtertype_;
            INTERVALFILTER intervalfilter_;
            //int intervalfilter_;
            tgt::vec3   pos_;
            tgt::vec3   normal_;
            float       radius_;
            tgt::Bounds bounds_;
            tgt::vec2   interval_;
			bool allDarkMatter_;
			bool allBaryons_;
			bool starBaryons_;
			bool windBaryons_;
			bool gasBaryons_;
			bool agn_;

            
      /**
       * "True" variables for converting the relative values contained in the filter into
       * exact values based on the current timestep. Necessary for handling comoving 
       * coordinates.
       */
       
            tgt::vec3   truepos_;
            float       trueradius_;
            tgt::Bounds truebounds_;
            tgt::vec2   trueinterval_;            
            
      /**
       * Custom serialization and deserialization to include the PartFilter objects.
       */

            virtual void serialize(Serializer& s) const override;
            virtual void deserialize(Deserializer& s) override;

        };
        std::vector<PartFilter> filters_;


        // Ports
        CMParticlePort                  inport_;
        CMParticlePort                  outport_;

        // Current timestep
        FloatProperty                   timeStep_;
        FloatProperty                   linkedTimeStep_;
		BoolProperty                    useLinkedTime_;

        // Properties for filter control GUI
        BoolProperty                    enabled_;
		BoolProperty                    anyTimeFilter_;
        BoolProperty                    invert_;
        IntProperty                     currentfilter_;
        OptionProperty<FILTERTYPE>      filtertype_;
        OptionProperty<FILTERCONNECTION> filterconnection_;
        OptionProperty<INTERVALFILTER>  intervalfilter_;
        FloatVec3Property               pos_;
        FloatVec3Property               normal_;
        FloatProperty                   radius_;
        FloatBoundingBoxProperty        bounds_;
        FloatVec2Property               interval_;
		//PropertyVector					typeSelection_;
		BoolProperty                   allDarkMatter_;
		BoolProperty                   allBaryons_;
		BoolProperty                   starBaryons_;
		BoolProperty                   windBaryons_;
		BoolProperty                   gasBaryons_;
		BoolProperty                   agn_;


        BoolProperty                    autoEval_;
        ButtonProperty                  evalButton_;
        ButtonProperty                  addFilter_;
        ButtonProperty                  removeFilter_;
        


		// Adds a new empty filter at the end of the vector filters_ 
        void addFilter();
        
        // Deletes the last filter in filters_
        void removeLastFilter();
        
        // Copies the values of the currently selected vector element into the "active filter"
        void popFilter();
        
        // Copies the "active filter" values into the currently selected vector element
        void pushFilter();

		//filter at timestep
		void doFilter(CMParticleDataTimeSlice* particleData, char* enabledState, int k, bool connectByAnd);
        
        const CMParticleData* particleDataIn_;
        CMParticleData* particleData_;
        std::vector<CMParticleDataTimeSlice*> timeSlices_;
        const CMParticleDataTimeSlice* timeSlicePointer_;
        //CMParticleDataTimeSlice timeSlice_;
        std::vector<CMParticle> particles_;

		// Adds the current filter to the vector and activates the filters for use in process()
        void onEvalButtonPressed();
		
		// Handles access of uninitialized/illegal positions in the filters_ vector
        void switchCurrentFilter();
        
        // Handles UI element visibility based on type of current filter
        void switchFilterType();

		// Handles checked boxes on type filter
		void onTypeChanged(TYPECHANGED box);

        bool doFilterOnNextEvaluation_;

		void blockTypeCallBacks(bool block) {
			allDarkMatter_.blockCallbacks(block);
			allBaryons_.blockCallbacks(block);
			starBaryons_.blockCallbacks(block);
			windBaryons_.blockCallbacks(block);
			gasBaryons_.blockCallbacks(block);
			agn_.blockCallbacks(block);
		}
        

    };
}   // end of namespace
#endif