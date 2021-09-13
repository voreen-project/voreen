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
#include "tgt/filesystem.h"
#include "cosmologydatasource.h"
#include "md5/md5.h"
#include "voreen/core/voreenapplication.h"

namespace voreen{

namespace{
class MemoryMapTimeSlice : public CMParticleDataTimeSlice{
public:
    MemoryMapTimeSlice(tgt::Bounds universeBounds);
    virtual ~MemoryMapTimeSlice();

    virtual const CMParticle* startUsingParticles() override;
    virtual void              finishedUsingParticles() override;
    virtual int               getNumberOfParticles() override;
    virtual float             getTimeStep() override;
    virtual const int*        getRemappingTable() override;
    virtual tgt::Bounds       getBounds() override;
    virtual float             geth0() override;


    void                      loadRemappingTable();
    std::string               getCacheEntryName(std::string nameString);

    boost::iostreams::mapped_file_source file_;
    int                numberOfParticles_;
    float              timeStep_;
    std::string        filename_;

    int* remapping_;
    static const std::string loggerCat_; ///< category used in logging
};
const std::string MemoryMapTimeSlice::loggerCat_("voreen.MemoryMappedTimeSlice");

MemoryMapTimeSlice::MemoryMapTimeSlice(tgt::Bounds universeBounds)
    : CMParticleDataTimeSlice(universeBounds)
{
    remapping_ = 0;
}

MemoryMapTimeSlice::~MemoryMapTimeSlice(){
    delete[] remapping_;
}

const CMParticle* MemoryMapTimeSlice::startUsingParticles(){
    if (!file_.is_open())
        file_.open(filename_);
    const char* d = (file_.data());
    return reinterpret_cast<const CMParticle*>(d);
}

void MemoryMapTimeSlice::finishedUsingParticles(){
    file_.close();
}

int MemoryMapTimeSlice::getNumberOfParticles(){
    return numberOfParticles_;    
}

float MemoryMapTimeSlice::getTimeStep(){
    return timeStep_;
}

const int* MemoryMapTimeSlice::getRemappingTable(){
    if (!remapping_){
        //loadRemappingTable();
		remapping_ = new int[numberOfParticles_];

		for (int i = 0; i < numberOfParticles_; i++) {
			remapping_[i] = i;
		}
    }

	return remapping_;
    //return 0;
}

tgt::Bounds MemoryMapTimeSlice::getBounds() {
    tgt::Bounds b = getUniverseBounds();
    tgt::mat4 normMatrix = getNormalisationTransformation();
    tgt::vec4 llf1(b.getLLF(), 1);
    tgt::vec4 urb1(b.getURB(), 1);
    return tgt::Bounds(normMatrix*llf1.xyz(), normMatrix*urb1.xyz());


}
float MemoryMapTimeSlice::geth0() {
    return 1.f;
}

void MemoryMapTimeSlice::loadRemappingTable(){
    std::string cacheFileName = getCacheEntryName(filename_);
    const std::string cacheFilePath = VoreenApplication::app()->getCachePath(cacheFileName);
    const std::string cacheFolderPath = VoreenApplication::app()->getCachePath("");
    if(!tgt::FileSystem::dirExists(cacheFolderPath)) {
        if(!tgt::FileSystem::createDirectoryRecursive(cacheFolderPath)) {
            //LERROR("Error encountered while loading halo cache file: Could not create voreen cache directory.");
            //mergerTree_ = nullptr;
            return;
        }
    }
    if(!tgt::FileSystem::fileExists(cacheFilePath)) {
        LDEBUG("Store remapping table " + cacheFileName);
        remapping_ = buildRemappingTable();
        std::ofstream stream(cacheFilePath, std::ofstream::binary);
        stream.write((char*)remapping_, sizeof(int)*getNumberOfParticles());
        return;
    }else{
        LDEBUG("Loaded remapping table" + cacheFileName);
        std::ifstream stream(cacheFilePath, std::ofstream::binary);
        remapping_ = new int[getNumberOfParticles()];
        stream.read((char*)remapping_, sizeof(int)*getNumberOfParticles());
        return;
    }
}

std::string MemoryMapTimeSlice::getCacheEntryName(std::string nameString){
    MD5_CTX md5ctx;
    MD5_Init(&md5ctx);
    MD5_Update(&md5ctx, (void*)nameString.c_str(), nameString.size());
    unsigned char md5[16];
    char md5str[33];
    MD5_Final(md5, &md5ctx);

    //MD5_Final(md5, &md5ctx);
    int i=0;
    for(i=0; i<16; ++i) {
        sprintf(md5str+i*2, "%02x", md5[i]);
    }
    md5str[32] = 0;
    LDEBUG("Hash for " + nameString + " is " + md5str);
    return std::string(md5str);
}

}

const std::string CosmologyDataSource::loggerCat_("voreen.CosmologyDataSource");

CosmologyDataSource::CosmologyDataSource()
    : Processor()
    , outport_(Port::OUTPORT, "particlehandle.output", "Particle Data Output")
    , fileFolder_("particleURL", "Load Particles", "Folder for Particle Files", "", "", FileDialogProperty::DIRECTORY)
{
    addPort(outport_);
    addProperty(fileFolder_);
}

CosmologyDataSource::~CosmologyDataSource(){
}

Processor* CosmologyDataSource::create() const{
    return new CosmologyDataSource;
}

void CosmologyDataSource::initialize() {
}

void CosmologyDataSource::deinitialize() {
    outport_.setData(nullptr);
}

void CosmologyDataSource::process(){
    CMParticleData* particleData = readParticleData();
    if (particleData){
        outport_.setData(particleData);
    }
}

CMParticleDataTimeSlice* CosmologyDataSource::tryReadFile(const std::string& filestring) {
    MemoryMapTimeSlice * timeSlice = new MemoryMapTimeSlice(tgt::Bounds()); //dummy bounds, will be set later
    try{
        //timeSlice->file_.open(filestring);

		size_t numberoffset = filestring.find("Full.cosmo.f");
		std::string number = std::string(filestring).erase(0, numberoffset + 12);
		int timestep = std::stoi(number);
        
		// hardcoded as size does not change
		timeSlice->numberOfParticles_ = 523519;
		//timeSlice->numberOfParticles_  = timeSlice->file_.size() / 64;
        timeSlice->timeStep_           = (float)timestep;
        //R0 read from file is in comoving coordinates => multiply by timeStep_ to get actual boxradius
        float boxRadius                = 32.f;
        timeSlice->universeBounds_     = tgt::Bounds(tgt::vec3(-boxRadius), tgt::vec3(boxRadius));
        timeSlice->filename_           = filestring;

		

        //timeSlice->finishedUsingParticles();
    }catch(...){
        delete timeSlice;
        timeSlice = 0;
        throw;
    }

    return timeSlice;
}

bool compareTimeStep(CMParticleDataTimeSlice* t1, CMParticleDataTimeSlice* t2)
{
	return (t1->getTimeStep() < t2->getTimeStep());
}

CMParticleData* CosmologyDataSource::readParticleData() {
    std::vector<CMParticleDataTimeSlice*> timeSlices;

    const std::string directoryString = fileFolder_.get();
    std::vector<std::string> filesInDir = tgt::FileSystem::readDirectory(directoryString, false, false);
    if (filesInDir.size() == 0){
        return nullptr;
    }

    for(auto filename : filesInDir) {
        try {
            timeSlices.push_back(tryReadFile(directoryString + "/" + filename));
        } catch (std::invalid_argument& e)
        {
            LWARNING("Error while reading " + filename + ": " + e.what());
        } catch(boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<std::ios_base::failure> > e) {
            LWARNING("Could not open file " + filename + ": " + e.what());
        }
    }
    //We could not find any SDFiles
    if (timeSlices.empty()){
        return nullptr;
	}
	else {
		std::sort(timeSlices.begin(), timeSlices.end(), compareTimeStep);
	}

    return new CMParticleData(timeSlices);
}
}
