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
#include "cmhalodatasource.h"
#include "../datastructures/cmmergertree.h"
#include "../datastructures/cmparticledata.h"
#include "../datastructures/cmmergertree.h"
#include "../datastructures/cmparticledata.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/stringutils.h"
#include "tgt/filesystem.h"
#include "md5/md5.h"

#ifdef HDF5_FOUND
#include "H5Cpp.h"
#endif




namespace voreen{


static void ignoreProperty(std::stringstream& sstream, int n) {
    float trash;
    for(int i=0; i<n; i++){
        sstream>> trash;
    }

}

CMHaloDataSource::CMHaloDataSource()
    : Processor()
      , outport_(Port::OUTPORT, "halohandle.output", "Halo Data Output")
      , fileFolder_("haloURL", "Load Halos", "Folder for Halo Files", "", "", FileDialogProperty::DIRECTORY)
      , galacticusInputFile_("galacticusInputFile", "Input File for Galacticus", "Input File for Galacticus", "", "*.xml", FileDialogProperty::SAVE_FILE)
      , galacticusOutputFile_("galacticusOutputFile", "Output File of Galacticus", "Output File of Galacticus", "", "*.hdf5", FileDialogProperty::OPEN_FILE)
      , galacticusExportButton_("galacticusExportButton_", "Export Galacticus Input File")
      , galacticusImportButton_("galacticusImportButton_", "Import Galacticus Output File")
      , mergerTreeFile_(nullptr)
      , mergerTree_(nullptr)
      , haloFilesLoaded_(false)
      , haloListsContainGalacticusData_(false)
      , treeValid_(true)
{
    addPort(outport_);
    addProperty(fileFolder_);
    addProperty(galacticusInputFile_);
    addProperty(galacticusExportButton_);
    addProperty(galacticusOutputFile_);
    addProperty(galacticusImportButton_);

    ON_CHANGE(fileFolder_, CMHaloDataSource, fileFolderChanged);
    ON_CHANGE(galacticusExportButton_, CMHaloDataSource, exportToGalacticus);
    ON_CHANGE(galacticusImportButton_, CMHaloDataSource, importFromGalacticus);
}

CMHaloDataSource::~CMHaloDataSource(){
}

Processor* CMHaloDataSource::create() const{
    return new CMHaloDataSource;
}

void CMHaloDataSource::initialize() {
}

void CMHaloDataSource::deinitialize() {
    closeFile();
    outport_.setData(nullptr, false);
}

void CMHaloDataSource::process(){
    outport_.setData(nullptr, false);
    closeFile();
    readmergerTreeData();
    outport_.setData(mergerTree_, false);
}

void CMHaloDataSource::closeFile() {
    if(mergerTreeFile_) {
        mergerTreeFile_->close();
        delete mergerTreeFile_;
        mergerTreeFile_ = nullptr;
    }
}
void CMHaloDataSource::fileFolderChanged() {
    const std::string directoryString = fileFolder_.get();
    if(!tgt::FileSystem::dirExists(directoryString)) {
        return;
    }
    //Calculate cache file name from md5-hash of source file names
    std::stringstream names;
    std::vector<std::string> filesInDir = tgt::FileSystem::readDirectory(directoryString, false, false);
    for(auto filename : filesInDir) {
        names << filename;
    }
    std::string nameString;
    names >> nameString;
    MD5_CTX md5ctx;
    MD5_Init(&md5ctx);
    MD5_Update(&md5ctx, (void*)nameString.c_str(), nameString.size());
    unsigned char md5[16];
    char md5str[33];
    MD5_Final(md5, &md5ctx);
    int i=0;
    for(i=0; i<16; ++i) {
        sprintf(md5str+i*2, "%02x", md5[i]);
    }
    md5str[32] = 0;
    haloFilesIdentifier_ = std::string(md5str);

    //Get rid of old halo lists
    haloLists_.clear();
    haloFilesLoaded_ = false;
    haloListsContainGalacticusData_ = false;
}
void CMHaloDataSource::invalidateMergerTree() {
    treeValid_ = false;
}
void CMHaloDataSource::exportToGalacticus() {
    const std::string directoryString = fileFolder_.get();
    const std::string hdfFile = tgt::FileSystem::dirName(galacticusInputFile_.get())
        + "/voreenOut" + haloFilesIdentifier_ + ".hdf5";
    std::cout << hdfFile << std::endl;
    std::vector<std::vector<CMHalo>>& haloLists = getHaloLists();
    writeGalacticusConfig(galacticusInputFile_.get(), hdfFile);
    exportToHDF5(hdfFile, haloLists);
}
void CMHaloDataSource::importFromGalacticus() {
    std::vector<std::vector<CMHalo>>& haloLists = getHaloLists();
    std::string galacticusFile = galacticusOutputFile_.get();
    if(!tgt::FileSystem::fileExists(galacticusFile)) {
        LERROR("Error encountered while loading halo cache file: Could not open halo cache file for reading");
        return;
    }
    if(!importFromHDF5(galacticusFile, haloLists)) {
        return;
    }
    haloListsContainGalacticusData_ = true;
    //We added new data fields to the halo files and thus have to update the tree
    invalidateMergerTree();
}

void CMHaloDataSource::readmergerTreeData() {
    std::vector<CMTimeStepHaloList> stepList;
    const std::string directoryString = fileFolder_.get();
    const std::string cacheFilePath = getCacheFilePath();

    if(!tgt::FileSystem::dirExists(directoryString)) {
        mergerTree_ = nullptr;
        return;
    }
    const std::string cacheFolderPath = VoreenApplication::app()->getCachePath("");
    if(!tgt::FileSystem::dirExists(cacheFolderPath)) {
        if(!tgt::FileSystem::createDirectoryRecursive(cacheFolderPath)) {
            LERROR("Error encountered while loading halo cache file: Could not create voreen cache directory.");
            mergerTree_ = nullptr;
            return;
        }
    }
    if(!tgt::FileSystem::fileExists(cacheFilePath) || !treeValid_) {
        saveTreeToFile();
    }
    if(!tgt::FileSystem::fileExists(cacheFilePath)) {
        LERROR("Error encountered while loading halo cache file: Could not open halo cache file for reading");
        mergerTree_ = nullptr;
        return;
    }
    try {
        mergerTreeFile_ = new boost::iostreams::mapped_file_source(cacheFilePath);
    } catch(...) {
        LERROR("Error encountered while loading halo cache file: Could not open halo cache file for reading");
        mergerTree_ = nullptr;
        return;
    }
    mergerTree_ = reinterpret_cast<const CMMergerTree*>(mergerTreeFile_->data());
}
std::string CMHaloDataSource::getCacheFilePath() const {
    return VoreenApplication::app()->getCachePath(haloFilesIdentifier_);
}
void CMHaloDataSource::saveTreeToFile() {
    std::vector<std::vector<CMHalo>>& hlists = getHaloLists();
    int numberOfHalos = 0;
    for(const auto& hlist : hlists) {
        numberOfHalos += hlist.size();
    }

    const size_t treeBlockSize = sizeof(CMMergerTree);
    const size_t timestepBlockSize = hlists.size()*sizeof(CMTimeStepHaloList);
    const size_t haloBlockSize = numberOfHalos*sizeof(CMHalo);
    boost::iostreams::mapped_file_params fileParams;
    fileParams.path = getCacheFilePath();
    fileParams.new_file_size = treeBlockSize+timestepBlockSize+haloBlockSize;
    try {
        boost::iostreams::mapped_file_sink outfile(fileParams);
        CMMergerTree* tree = reinterpret_cast<CMMergerTree*>(outfile.data());
        CMTimeStepHaloList* timestep = reinterpret_cast<CMTimeStepHaloList*>(outfile.data()+treeBlockSize);
        CMHalo* halo = reinterpret_cast<CMHalo*>(outfile.data()+treeBlockSize+timestepBlockSize);
#define BYTE_DIST(a, b) std::distance(reinterpret_cast<char*>(a), reinterpret_cast<char*>(b))
        tree->halosBeginOffset = BYTE_DIST(tree, halo);
        tree->hlistBeginOffset = BYTE_DIST(tree, timestep);
        tree->containsGalacticusData_ = haloListsContainGalacticusData_;
        for(auto& hlist : hlists){
            timestep->beginOffset = BYTE_DIST(timestep, halo);
            for(auto& h: hlist){
                *halo = h;
                ++halo;
            }
            //Just take the time of last halo in list as they should all have the same value for a
            timestep->a = (halo-1)->scale;
            timestep->endOffset = BYTE_DIST(timestep, halo);
            ++timestep;
        }
        tree->hlistEndOffset = BYTE_DIST(tree, timestep);
        tree->halosEndOffset = BYTE_DIST(tree, halo);
#undef BYTE_DIST
    } catch(...) {
        LERROR("Error encountered while writing halo cache file.");
        return;
    }
    treeValid_ = true;
}
void CMHaloDataSource::readHaloListFile(const std::string& fileName, std::vector<CMHalo>& halos) const {
    tgt::RegularFile file(fileName);

    while(!file.eof()) {
        std::string line=file.getLine();
        std::stringstream sstream;

        if(line[0]!='#') {
            CMHalo ha;

            sstream << line;

            sstream >> ha.scale;
            sstream >>ha.ID;
            ha.origID = ha.ID;
            tgtAssert(ha.ID >= 0, "id is negative");
            sstream >>ha.descendantScale;
            sstream >>ha.descendantID;
            ignoreProperty(sstream,1);
            sstream >>ha.hostID;
            sstream >>ha.rootHostID;
            ignoreProperty(sstream,3);
            sstream >>ha.mass;
            sstream >>ha.radius;
            sstream >>ha.scaleRadius;
            ignoreProperty(sstream,4);
            sstream >>ha.pos.x;
            sstream >>ha.pos.y;
            sstream >>ha.pos.z;
            sstream >>ha.velocity.x;
            sstream >>ha.velocity.y;
            sstream >>ha.velocity.z;
            sstream >>ha.angularMomenta.x;
            sstream >>ha.angularMomenta.y;
            sstream >>ha.angularMomenta.z;
            sstream >>ha.spinParameter;
            ignoreProperty(sstream,2);
            sstream >>ha.origTreeRootID;
            ha.parentID = CMMergerTree::NO_HALO_ID;
            ha.spouseID = CMMergerTree::NO_HALO_ID;
            ha.satelliteID = CMMergerTree::NO_HALO_ID;
            ha.siblingSatelliteID = CMMergerTree::NO_HALO_ID;
            halos.push_back(ha);
        }
    }
}
static void calcTreeLinks(std::vector<std::vector<CMHalo>>& hlists) {
    std::vector<CMHalo*> haloByOldID;
    std::vector<CMHalo*> haloByNewID;
    int maxID = 0;
    int idcounter = 0;
    //Find out maximum of old ids which might be spread out
    for(auto& hlist : hlists){
        for(auto& halo : hlist){
            if(maxID<halo.ID){
                maxID = halo.ID;
            }
        }
    }
    //Set new id and build up LUT newIDByOldID
    haloByOldID.resize(maxID+1);
    haloByNewID.resize(maxID+1);
    for(auto& hlist : hlists){
        for(auto& halo : hlist){
            haloByOldID.at(halo.ID) = &halo;
            halo.ID = idcounter++;
            haloByNewID.at(halo.ID) = &halo;
        }
    }
    //Convert existing old ids to new ids
    for(auto& hlist : hlists){
        for(auto& halo : hlist){
            if(halo.descendantID!=CMMergerTree::NO_HALO_ID) {
                halo.descendantID = haloByOldID[halo.descendantID]->ID;
            }
            if(halo.hostID!=CMMergerTree::NO_HALO_ID) {
                halo.hostID = haloByOldID[halo.hostID]->ID;
            }
            if(halo.rootHostID!=CMMergerTree::NO_HALO_ID) {
                halo.rootHostID = haloByOldID[halo.rootHostID]->ID;
            }
        }
    }
    //Calculate parentID and spouseID backreferences
    for(int i=0; i<static_cast<int>(hlists.size())-1;i++){
        for(auto& h: hlists[i]){
            CMHalo* descHalo = haloByNewID[h.descendantID];
            h.spouseID = descHalo->parentID;
            descHalo->parentID= h.ID;
        }
    }
    //Calculate satelliteID and nextSatelliteID backreferences
    for(auto& hlist : hlists){
        for(auto& h: hlist){
            if(h.hostID != CMMergerTree::NO_HALO_ID) {
                CMHalo* hostHalo = haloByNewID[h.hostID];
                h.siblingSatelliteID = hostHalo->satelliteID;
                hostHalo->satelliteID = h.ID;
            }
        }
    }
}

std::vector<std::vector<CMHalo>>& CMHaloDataSource::getHaloLists() {
    if(!haloFilesLoaded_) {
        //std::cout << "Hlist get start" << std::endl;
        const std::string directoryString = fileFolder_.get();
        std::vector<std::string> filesInDir = tgt::FileSystem::readDirectory(directoryString, false, false);
        for(auto filename : filesInDir) {
            haloLists_.emplace_back();
            readHaloListFile(directoryString + "/" + filename, haloLists_.back());
        }
        std::sort(haloLists_.begin(),haloLists_.end(),[&](std::vector<CMHalo> h1, std::vector<CMHalo> h2){return h1.at(0).scale < h2.at(0).scale;});
        calcTreeLinks(haloLists_);
        haloFilesLoaded_ = true;
        //std::cout << "Hlist get end" << std::endl;
    }
    return haloLists_;
}
#ifdef HDF5_FOUND
void createAndFillAttributeF64(const H5::H5Location& fg, const char* name, const double data) {
    static const H5::DataSpace dataspace; //scalar data space
    H5::Attribute attribute = fg.createAttribute(name, H5::PredType::IEEE_F64LE, dataspace);
    attribute.write(H5::PredType::IEEE_F64LE, &data);
}
void createAndFillAttributeI64(const H5::H5Location& fg, const char* name, const int64_t data) {
    static const H5::DataSpace dataspace; //scalar data space
    H5::Attribute attribute = fg.createAttribute(name, H5::PredType::STD_I64LE, dataspace);
    attribute.write(H5::PredType::STD_I64LE, &data);
}
void createAndFillAttributeI32(const H5::H5Location& fg, const char* name, const int32_t data) {
    static const H5::DataSpace dataspace; //scalar data space
    H5::Attribute attribute = fg.createAttribute(name, H5::PredType::STD_I32LE, dataspace);
    attribute.write(H5::PredType::STD_I32LE, &data);
}
void createAndWriteToDataSet(const H5::CommonFG& fg, const char* name, const H5::DataType& type, const H5::DataSpace& dataspace, const void* data) {
    H5::DataSet dataset = fg.createDataSet(name, type, dataspace);
    dataset.write(data, type);
}
void copyValue(int64_t* arr, const int64_t i, const int val) {
    arr[i] = val;
}
void copyValue(double* arr, const int64_t i, const float val) {
    arr[i] = val;
}
void copyValue(double* arr, const int64_t i, const tgt::vec3& vec) {
    arr[3*i  ] = vec.x;
    arr[3*i+1] = vec.y;
    arr[3*i+2] = vec.z;
}
class TreeRootsIndex {
    public:
    TreeRootsIndex(const std::vector<CMHalo>& treeRoots)
        : treeRoots_(treeRoots)
    {
        auto minmaxID = std::minmax_element(treeRoots_.begin(), treeRoots_.end(), [] (CMHalo h1, CMHalo h2) { return h1.origID < h2.origID; });
        smallestID_ = minmaxID.first->origID;
        largestID_ = minmaxID.second->origID;
        rootByID_ = new int[largestID_-smallestID_+1];
        for(int i=0; i<treeRoots_.size(); ++i) {
            rootByID_[treeRoots[i].origID-smallestID_] = i;
        }
    }
    ~TreeRootsIndex() {
        delete[] rootByID_;
    }
    int treeIDByRootID(int i) {
        if(i < smallestID_ or i > largestID_) {
            LERROR("TreeRootsIndex.treeIDByRootID: Invalid parameter id: " << smallestID_ << " < " << i << " < " << largestID_);
        }
        return rootByID_[i-smallestID_];
    }

    private:
    const std::vector<CMHalo>& treeRoots_;
    int smallestID_;
    int largestID_;
    int* rootByID_;
    static const std::string loggerCat_; ///< category used in logging
};
const std::string TreeRootsIndex::loggerCat_("voreen.treerootsindex") ;
class H5MergerTree {
public:
    std::vector<int64_t> IDs_;
    std::vector<int64_t> hostIDs_;
    std::vector<int64_t> descendentIDs_;
    std::vector<double> scaleRadii_;
    std::vector<double> masses_;
    std::vector<double> positions_;
    std::vector<double> spins_;
    std::vector<double> times_;
    std::vector<double> velocities_;
    std::vector<double> angularMomenta_;

    void add(const CMHalo& h) {
        IDs_.push_back(h.ID);
        hostIDs_.push_back((h.hostID == -1) ? h.ID : h.hostID);
        descendentIDs_.push_back(h.descendantID);
        scaleRadii_.push_back(h.scaleRadius/1000); //Radii are given in kpc/h, but we want Mpc/h
        masses_.push_back(h.mass);
        positions_.push_back(h.pos.x);
        positions_.push_back(h.pos.y);
        positions_.push_back(h.pos.z);
        spins_.push_back(h.spinParameter);
        times_.push_back(h.scale);
        velocities_.push_back(h.velocity.x);
        velocities_.push_back(h.velocity.y);
        velocities_.push_back(h.velocity.z);
        angularMomenta_.push_back(h.angularMomenta.x);
        angularMomenta_.push_back(h.angularMomenta.y);
        angularMomenta_.push_back(h.angularMomenta.z);
    }
    void add(const H5MergerTree& t) {
        IDs_.insert(IDs_.end(), t.IDs_.begin(), t.IDs_.end());
        hostIDs_.insert(hostIDs_.end(), t.hostIDs_.begin(), t.hostIDs_.end());
        descendentIDs_.insert(descendentIDs_.end(), t.descendentIDs_.begin(), t.descendentIDs_.end());
        scaleRadii_.insert(scaleRadii_.end(), t.scaleRadii_.begin(), t.scaleRadii_.end());
        masses_.insert(masses_.end(), t.masses_.begin(), t.masses_.end());
        positions_.insert(positions_.end(), t.positions_.begin(), t.positions_.end());
        spins_.insert(spins_.end(), t.spins_.begin(), t.spins_.end());
        times_.insert(times_.end(), t.times_.begin(), t.times_.end());
        velocities_.insert(velocities_.end(), t.velocities_.begin(), t.velocities_.end());
        angularMomenta_.insert(angularMomenta_.end(), t.angularMomenta_.begin(), t.angularMomenta_.end());
    }
    hsize_t size() {
        return IDs_.size();
    }
};
#endif
void CMHaloDataSource::exportToHDF5(const std::string& exportFilePath, std::vector<std::vector<CMHalo>>& haloLists) const {
#ifdef HDF5_FOUND
    std::sort(haloLists.begin(),haloLists.end(),[](std::vector<CMHalo> h1, std::vector<CMHalo> h2){return h1.at(0).scale < h2.at(0).scale;});
    //std::random_shuffle(haloLists.begin(),haloLists.end());
    hsize_t numberOfHalos = 0;
    hsize_t numberOfTrees = 0;
    const std::vector<CMHalo>* treeRoots = nullptr;
    for(const std::vector<CMHalo>& halos : haloLists) {
        numberOfHalos += halos.size();
        if(halos[0].scale == 1.0f) {
            treeRoots = &halos;
        }
    }
    assert(treeRoots != nullptr);
    numberOfTrees = treeRoots->size();
    // Create the data space for the dataset.
    H5::DataSpace scalar_dataspace(1, &numberOfHalos);
    //H5::DataSpace tree_dataspace(1, &numberOfTrees);
    hsize_t one = 1;
    H5::DataSpace tree_dataspace(1, &one);
    hsize_t dims[2];
    dims[0] = numberOfHalos;
    dims[1] = 3;
    H5::DataSpace vec3_dataspace(2, dims);

    int32_t* firstNode     = new int32_t[numberOfTrees];
    int32_t* numberOfNodes = new int32_t[numberOfTrees];
    int64_t* treeIndex     = new int64_t[numberOfTrees];

    TreeRootsIndex treeRootsIndex(*treeRoots);
    std::vector<H5MergerTree> trees(numberOfTrees);

    for(const std::vector<CMHalo>& halos : haloLists) {
        for(const CMHalo& h : halos) {
            int treeID = treeRootsIndex.treeIDByRootID(h.origTreeRootID);
            trees.at(treeID).add(h);
        }
    }
    H5MergerTree allTrees;
    int i=0;
    for(H5MergerTree& t : trees) {
        allTrees.add(t);
        treeIndex[i] = i;
        firstNode[i] = allTrees.size()-t.size();;
        numberOfNodes[i] = t.size();
        ++i;
    }
    firstNode[0] = 0;
    numberOfNodes[0] = numberOfHalos;
    treeIndex[0] = 0;
    // Try block to detect exceptions raised by any of the calls inside it
    try
    {
        // Turn off the auto-printing when failure occurs so that we can
        // handle the errors appropriately
        //H5::Exception::dontPrint();

        // Create a new file using the default property lists.
        H5::H5File file(exportFilePath, H5F_ACC_TRUNC);


        // Create Groups and write datasets
        {
            H5::Group group = file.createGroup("cosmology");

            //createAndFillAttributeF64(group, "HubbleParam", 0.68); //woher?
        }
        //Group finder?
        {
            H5::Group group = file.createGroup("haloTrees");

            //?
            createAndFillAttributeI32(group, "haloMassesIncludeSubhalos", 0);
            createAndFillAttributeI32(group, "haloAngularMomentaIncludeSubhalos", 0);
            //createAndFillAttributeI32(group, "treesAreSelfContained", 0);
            //createAndFillAttributeI32(group, "treesHaveSubhalos", 0);
            createAndFillAttributeI32(group, "velocitiesIncludeHubbleFlow", 0);
            createAndFillAttributeI32(group, "positionsArePeriodic", 1);

            createAndWriteToDataSet(group, "descendentIndex", H5::PredType::STD_I64LE , scalar_dataspace, allTrees.descendentIDs_.data());
            createAndWriteToDataSet(group, "scaleRadius"    , H5::PredType::IEEE_F64LE, scalar_dataspace, allTrees.scaleRadii_.data()    );
            createAndWriteToDataSet(group, "hostIndex"      , H5::PredType::STD_I64LE , scalar_dataspace, allTrees.hostIDs_.data()       );
            createAndWriteToDataSet(group, "nodeIndex"      , H5::PredType::STD_I64LE , scalar_dataspace, allTrees.IDs_.data()           );
            createAndWriteToDataSet(group, "nodeMass"       , H5::PredType::IEEE_F64LE, scalar_dataspace, allTrees.masses_.data()        );
            //createAndWriteToDataSet(group, "particleCount"  , H5::PredType::STD_I64LE , scalar_dataspace, ??         );
            //createAndWriteToDataSet(group, "particleStart"  , H5::PredType::STD_I64LE , scalar_dataspace, ??         );
            createAndWriteToDataSet(group, "position"       , H5::PredType::IEEE_F64LE, vec3_dataspace  , allTrees.positions_.data()     );
            createAndWriteToDataSet(group, "spin"           , H5::PredType::IEEE_F64LE, scalar_dataspace, allTrees.spins_.data()         );
            //Maybe instead of "time": "redshift" or "expansionFactor"?
            createAndWriteToDataSet(group, "expansionFactor", H5::PredType::IEEE_F64LE, scalar_dataspace, allTrees.times_.data()         );
            createAndWriteToDataSet(group, "velocity"       , H5::PredType::IEEE_F64LE, vec3_dataspace  , allTrees.velocities_.data()    );
            createAndWriteToDataSet(group, "angularMomentum", H5::PredType::IEEE_F64LE, vec3_dataspace  , allTrees.angularMomenta_.data());
        }
        //Mergertrees? optional apparently
        //particles?
        //provenance?
        {
            H5::Group group = file.createGroup("simulation");
            createAndFillAttributeF64(group, "boxSize", 62.5);
        }
        {
            H5::Group group = file.createGroup("treeIndex");

            createAndWriteToDataSet(group, "firstNode"     , H5::PredType::STD_I32LE, tree_dataspace, firstNode    );
            createAndWriteToDataSet(group, "numberOfNodes" , H5::PredType::STD_I32LE, tree_dataspace, numberOfNodes);
            createAndWriteToDataSet(group, "treeIndex"     , H5::PredType::STD_I64LE, tree_dataspace, treeIndex    );
        }
        {
            H5::Group group = file.createGroup("units");

            //Only for our correct for the given IEEE SciVis Contest 2015 dataset:
            //Lengths are in Mpc/h (see Galacticus manual v0.9.3 page 2006
            createAndFillAttributeI32(group, "lengthHubbleExponent", -1);
            createAndFillAttributeI32(group, "lengthScaleFactorExponent", 1);
            createAndFillAttributeF64(group, "lengthUnitsInSI", 3.08568e+22);

            //Velocities are in km/s
            createAndFillAttributeI32(group, "velocityHubbleExponent", 0);
            createAndFillAttributeI32(group, "velocityScaleFactorExponent", 0);
            createAndFillAttributeF64(group, "velocityUnitsInSI", 1000);

            //Masses are in Msun/h
            createAndFillAttributeI32(group, "massHubbleExponent", -1);
            createAndFillAttributeI32(group, "massScaleFactorExponent", 0);
            //M stands for Mega?
            //createAndFillAttributeF64(group, "massUnitsInSI", 1.98892e+36);
            //or: M stands for mass?
            createAndFillAttributeF64(group, "massUnitsInSI", 1.98892e+30);

            //Who knows?
            createAndFillAttributeI32(group, "timeHubbleExponent", 0);
            createAndFillAttributeI32(group, "timeScaleFactorExponent", 0);
            createAndFillAttributeF64(group, "timeUnitsInSI", 3.1556926e+16);
        }

    }  // end of try block

    // catch failure caused by the H5File operations
    catch(H5::FileIException error)
    {
        error.printError();
    }

    // catch failure caused by the DataSet operations
    catch(H5::DataSetIException error)
    {
        error.printError();
    }

    // catch failure caused by the DataSpace operations
    catch(H5::DataSpaceIException error)
    {
        error.printError();
    }

    delete[] firstNode;
    delete[] numberOfNodes;
    delete[] treeIndex;

#endif
}
#ifdef HDF5_FOUND
template<typename dataType>
void readDataset(dataType*& data, hsize_t& numElements, const H5::CommonFG& file, const char* path) {
    H5::DataType h5type;
    if(std::is_same<dataType, double>::value) {
        h5type = H5::PredType::IEEE_F64LE;
    } else if(std::is_same<dataType, int32_t>::value) {
        h5type = H5::PredType::STD_I32LE;
    } else if(std::is_same<dataType, int64_t>::value) {
        h5type = H5::PredType::STD_I64LE;
    }
    H5::DataSet dataset = file.openDataSet(path);
    numElements = dataset.getSpace().getSelectNpoints();
    data = new dataType[numElements];
    dataset.read(data, h5type);
}
#endif
bool CMHaloDataSource::importFromHDF5(const std::string& importFilePath, std::vector<std::vector<CMHalo>>& hlists) const {
#ifdef HDF5_FOUND

    int maxID = 0;
    //Find out maximum of old ids which might be spread out
    for(auto& hlist : hlists){
        for(auto& halo : hlist){
            if(maxID<halo.origID){
                maxID = halo.ID;
            }
        }
    }
    //Build up LUT haloByID
    std::vector<CMHalo*> haloByID(maxID+1, nullptr);
    for(auto& hlist : hlists){
        for(auto& halo : hlist){
            haloByID.at(halo.ID) = &halo;
        }
    }

    try
    {
        int64_t* ID = nullptr;
        double* blackHoleMass;
        double* blackHoleSpin;
        double* spheroidRadius;
        double* spheroidMassGas;
        double* spheroidVelocity;
        double* diskRadius;
        double* diskMassGas;
        double* diskVelocity;

        hsize_t numberOfNodes;
        hsize_t lastNumElements;
        // Turn off the auto-printing when failure occurs so that we can
        // handle the errors appropriately
        //H5::Exception::dontPrint();

        // Create a new file using the default property lists.
        H5::H5File file(importFilePath, H5F_ACC_RDONLY);
        H5::Group outputs = file.openGroup("/Outputs");
        for(int oid=1; oid <= outputs.getNumObjs(); ++oid) {
            //Relies on the specific galacticus output format, c++ h5 iterators are not implemented, yet
            std::string outputName = "Output" + std::to_string(oid);
            H5::Group currentOutput = outputs.openGroup(outputName);
            readDataset(ID              , numberOfNodes  , currentOutput, "nodeData/nodeIndex"        );
            readDataset(blackHoleMass   , lastNumElements, currentOutput, "nodeData/blackHoleMass"    );
            tgtAssert(numberOfNodes == lastNumElements, "Dataset size mismatch.");
            readDataset(blackHoleSpin   , lastNumElements, currentOutput, "nodeData/blackHoleSpin"    );
            tgtAssert(numberOfNodes == lastNumElements, "Dataset size mismatch.");
            readDataset(spheroidRadius  , lastNumElements, currentOutput, "nodeData/spheroidRadius"   );
            tgtAssert(numberOfNodes == lastNumElements, "Dataset size mismatch.");
            readDataset(spheroidMassGas , lastNumElements, currentOutput, "nodeData/spheroidMassGas"  );
            tgtAssert(numberOfNodes == lastNumElements, "Dataset size mismatch.");
            readDataset(spheroidVelocity, lastNumElements, currentOutput, "nodeData/spheroidVelocity" );
            tgtAssert(numberOfNodes == lastNumElements, "Dataset size mismatch.");
            readDataset(diskRadius      , lastNumElements, currentOutput, "nodeData/diskRadius"       );
            tgtAssert(numberOfNodes == lastNumElements, "Dataset size mismatch.");
            readDataset(diskMassGas     , lastNumElements, currentOutput, "nodeData/diskMassGas"      );
            tgtAssert(numberOfNodes == lastNumElements, "Dataset size mismatch.");
            readDataset(diskVelocity    , lastNumElements, currentOutput, "nodeData/diskVelocity"     );
            tgtAssert(numberOfNodes == lastNumElements, "Dataset size mismatch.");
            for(int i = 0; i < numberOfNodes; ++i) {
                CMHalo* halo = haloByID.at(ID[i]);
                tgtAssert(halo != nullptr, "Unknown halo ID.");
                halo->blackHoleMass = blackHoleMass[i];
                halo->blackHoleSpin = blackHoleSpin[i];
                halo->spheroidRadius = spheroidRadius[i];
                halo->spheroidMassGas = spheroidMassGas[i];
                halo->spheroidVelocity = spheroidVelocity[i];
                halo->diskRadius = diskRadius[i];
                halo->diskMassGas = diskMassGas[i];
                halo->diskVelocity = diskVelocity[i];
            }
            //std::cout << "Found " << numberOfNodes << " nodes." << std::endl;
            delete[] ID;
            delete[] blackHoleMass;
            delete[] blackHoleSpin;
            delete[] spheroidRadius;
            delete[] spheroidMassGas;
            delete[] spheroidVelocity;
            delete[] diskRadius;
            delete[] diskMassGas;
            delete[] diskVelocity;
        }
    }  // end of try block

    // catch failure caused by the H5File operations
    catch(H5::FileIException error)
    {
        error.printError();
        return false;
    }

    // catch failure caused by the DataSet operations
    catch(H5::DataSetIException error)
    {
        error.printError();
        return false;
    }

    // catch failure caused by the DataSpace operations
    catch(H5::DataSpaceIException error)
    {
        error.printError();
        return false;
    }
    catch(...) {
        LERROR("Unknown error while reading hdf5 file");
        return false;
    }
#endif
    return true;
}

std::string CMHaloDataSource::getRedshiftsString() {
    const std::vector<std::vector<CMHalo>>& haloLists = getHaloLists();
    std::ostringstream redshifts;
    redshifts << std::fixed;
    for(const auto& haloList : haloLists) {
        float expansionFactor = haloList.at(0).scale;
        float redshift = 1.0f/expansionFactor-1.0f;
        redshifts << redshift << " ";
    }
    return redshifts.str();
}
void CMHaloDataSource::writeGalacticusConfig(const std::string configFileName, const std::string& hdfFileName) {
    std::ofstream fs;
    fs.open(configFileName);
    fs << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
<parameters>\n\
    <formatVersion>2</formatVersion>\n\
    <version>0.9.4</version>\n\
    \n\
    <verbosityLevel value=\"3\"/>\n\
    \n\
    <!--  tree building options -->\n\
    <mergerTreeConstructMethod value=\"read\"/>\n\
    <mergerTreeReadFileName value=\"" + hdfFileName + "\"/>\n\
    <mergerTreeReadAllowSubhaloPromotions value=\"true\"/>\n\
    <mergerTreeReadPresetSpins value=\"true\"/>\n\
    <mergerTreeReadPresetUnphysicalSpins value=\"true\"/>\n\
    <treeNodeMethodPosition value=\"preset\"/>\n\
    <treeNodeMethodSatellite value=\"preset\"/>\n\
    <treeNodeMethodSpin value=\"preset\"/>\n\
    \n\
    <mergerTreeReadOutputTimeSnapTolerance value=\"1.0e-3\"/>\n\
    <outputRedshifts value=\""+ getRedshiftsString() +"\" />\n\
    \n\
    <!-- Component selection -->\n\
    <treeNodeMethodBasic value=\"standard\"/>\n\
    <treeNodeMethodBlackHole value=\"standard\"/>\n\
    <treeNodeMethodDarkMatterProfile value=\"scalePreset\"/>\n\
    <treeNodeMethodDisk value=\"exponential\"/>\n\
    <treeNodeMethodHotHalo value=\"standard\"/>\n\
    <treeNodeMethodSpheroid value=\"standard\"/>\n\
    <spheroidMassDistribution value=\"hernquist\"/>\n\
    \n\
    <!-- Cosmological parameters and options -->\n\
    <cosmologyFunctionsMethod value=\"matterLambda\"/>\n\
    <cosmologyParametersMethod value=\"simple\">\n\
    <HubbleConstant value=\"70.2\"/>\n\
    <OmegaMatter value=\"0.2725\"/>\n\
    <OmegaDarkEnergy value=\"0.7275\"/>\n\
    <OmegaBaryon value=\"0.0455\"/>\n\
    <temperatureCMB value=\"2.72548\"/>\n\
    </cosmologyParametersMethod>\n\
    <!-- Power spectrum options -->\n\
    <sigma_8 value=\"0.807\"/>\n\
    <transferFunctionMethod value=\"eisensteinHu1999\">  <neutrinoNumberEffective value=\"3.04\"/>\n\
    <neutrinoMassSummed value=\"0.0\"/>\n\
    </transferFunctionMethod>\n\
    <powerSpectrumMethod value=\"powerLaw\"/>\n\
    <powerSpectrumIndex value=\"0.961\"/>\n\
    <powerSpectrumReferenceWavenumber value=\"1.0\"/>\n\
    <powerSpectrumRunning value=\"0.0\"/>\n\
    <!-- Structure formation options -->\n\
    <linearGrowthMethod value=\"simple\"/>\n\
    <haloMassFunctionMethod value=\"Tinker2008\"/>\n\
    <criticalOverdensityMethod value=\"sphericalTopHat\"/>\n\
    <virialDensityContrastMethod value=\"sphericalCollapseMatterLambda\"/>\n\
    \n\
    <!-- Substructure hierarchy options -->\n\
    <nodeMergersMethod value=\"singleLevelHierarchy\"/>\n\
    \n\
    <!-- Dark matter halo structure options -->\n\
    <darkMatterProfileMethod value=\"NFW\"/>\n\
    <darkMatterProfileConcentrationMethod value=\"gao2008\"/>\n\
    <darkMatterProfileMinimumConcentration value=\"4\"/>\n\
    <haloSpinDistributionMethod value=\"Bett2007\"/>\n\
    <spinDistributionBett2007Alpha value=\"2.509\"/>\n\
    <spinDistributionBett2007Lambda0 value=\"0.04326\"/>\n\
    <randomSpinResetMassFactor value=\"2.0\"/>\n\
    \n\
    <!-- Halo accretion options -->\n\
    <accretionHalosMethod value=\"simple\"/>\n\
    <reionizationSuppressionRedshift value=\"10.5\"/>\n\
    <reionizationSuppressionVelocity value=\"35.0\"/>\n\
    \n\
    <!-- Hot halo ram pressure stripping options -->\n\
    <hotHaloRamPressureStrippingMethod value=\"Font2008\"/>\n\
    <hotHaloRamPressureForceMethod value=\"Font2008\"/>\n\
    <hotHaloRamPressureStrippingTimescaleMethod value=\"ramPressureAcceleration\"/>\n\
    <hotHaloOutflowStrippingEfficiency value=\"0.1\"/>\n\
    <hotHaloTrackStrippedGas value=\"true\"/>\n\
    \n\
    <!-- Galaxy ram pressure stripping options -->\n\
    <ramPressureStrippingMassLossRateDisksMethod value=\"null\"/>\n\
    <ramPressureStrippingMassLossRateSpheroidsMethod value=\"null\"/>\n\
    \n\
    <!-- Galaxy tidal stripping options -->\n\
    <tidalStrippingMassLossRateDisksMethod value=\"null\"/>\n\
    <tidalStrippingMassLossRateSpheroidsMethod value=\"null\"/>\n\
    <satellitesTidalFieldMethod value=\"null\"/>\n\
    \n\
    <!-- Galactic structure solver options -->\n\
    <galacticStructureRadiusSolverMethod value=\"adiabatic\"/>\n\
    <adiabaticContractionGnedinA value=\"0.73\"/>\n\
    <adiabaticContractionGnedinOmega value=\"0.7\"/>\n\
    <spheroidAngularMomentumAtScaleRadius value=\"0.5\"/>\n\
    \n\
    <!-- Galactic disk dynamics options -->\n\
    <barInstabilityMethod value=\"ELN\"/>\n\
    <stabilityThresholdGaseous value=\"0.7\"/>\n\
    <stabilityThresholdStellar value=\"1.1\"/>\n\
    \n\
    <!-- Star formation rate options -->\n\
    <starFormationTimescaleDisksMethod value=\"integratedSurfaceDensity\"/>\n\
    <starFormationRateSurfaceDensityDisksMethod value=\"KMT09\"/>\n\
    <molecularComplexClumpingFactorKMT09 value=\"5.0\"/>\n\
    <starFormationFrequencyKMT09 value=\"0.385\"/>\n\
    <molecularFractionFastKMT09 value=\"true\"/>\n\
    <starFormationDiskMinimumTimescale value=\"0.001\"/>\n\
    <starFormationTimescaleSpheroidsMethod value=\"dynamicalTime\"/>\n\
    <starFormationSpheroidEfficiency value=\"0.04\"/>\n\
    <starFormationSpheroidVelocityExponent value=\"2.0\"/>\n\
    <starFormationSpheroidMinimumTimescale value=\"0.001\"/>\n\
    \n\
    <!-- Stellar populations options -->\n\
    <stellarPopulationPropertiesMethod value=\"instantaneous\"/>\n\
    <stellarPopulationSpectraMethod value=\"Conroy-White-Gunn2009\"/>\n\
    <imfSelectionMethod value=\"fixed\"/>\n\
    <imfSelectionFixed value=\"Chabrier\"/>\n\
    <imfChabrierRecycledInstantaneous value=\"0.46\"/>\n\
    <imfChabrierYieldInstantaneous value=\"0.035\"/>\n\
    \n\
    <!-- AGN feedback options -->\n\
    <hotHaloExcessHeatDrivesOutflow value=\"true\"/>\n\
    <blackHoleHeatsHotHalo value=\"true\"/>\n\
    \n\
    <!-- Supernovae feedback options -->\n\
    <starFormationFeedbackDisksMethod value=\"powerLaw\"/>\n\
    <starFormationFeedbackSpheroidsMethod value=\"powerLaw\"/>\n\
    <diskOutflowVelocity value=\"250.0\"/>\n\
    <spheroidOutflowVelocity value=\"100.0\"/>\n\
    <diskOutflowExponent value=\"3.5\"/>\n\
    <spheroidOutflowExponent value=\"3.5\"/>\n\
    \n\
    <!-- Accretion disk properties -->\n\
    <accretionDisksMethod value=\"switched\"/>\n\
    <adafRadiativeEfficiencyType value=\"thinDisk\"/>\n\
    <accretionDiskSwitchedScaleAdafRadiativeEfficiency value=\"true\"/>\n\
    <accretionRateThinDiskMaximum value=\"0.30d0\"/>\n\
    <accretionRateThinDiskMinimum value=\"0.01d0\"/>\n\
    <adafAdiabaticIndex value=\"1.444\"/>\n\
    <adafEnergyOption value=\"pureADAF\"/>\n\
    <adafRadiativeEfficiency value=\"0.01\"/>\n\
    <adafViscosityOption value=\"fit\"/>\n\
    \n\
    <!-- Black hole options -->\n\
    <blackHoleBinaryMergersMethod value=\"Rezzolla2008\"/>\n\
    <blackHoleSeedMass value=\"100\"/>\n\
    <blackHoleWindEfficiency value=\"0.0024\"/>\n\
    <blackHoleWindEfficiencyScalesWithRadiativeEfficiency value=\"true\"/>\n\
    <bondiHoyleAccretionEnhancementHotHalo value=\"6.0\"/>\n\
    <bondiHoyleAccretionEnhancementSpheroid value=\"5.0\"/>\n\
    <bondiHoyleAccretionTemperatureSpheroid value=\"100\"/>\n\
    <bondiHoyleAccretionHotModeOnly value=\"true\"/>\n\
    \n\
    <!-- Satellite orbit options -->\n\
    <satelliteOrbitStoreOrbitalParameters value=\"true\"/>\n\
    \n\
    <!-- Galaxy merger options -->\n\
    <virialOrbitsMethod value=\"Benson2005\"/>\n\
    <satelliteMergingTimescalesMethod value=\"jiang2008\"/>\n\
    <mergingTimescaleMultiplier value=\"0.75\"/>\n\
    <satelliteMergingMassMovementsMethod value=\"simple\"/>\n\
    <minorMergerGasMovesTo value=\"spheroid\"/>\n\
    <satelliteMergingRemnantSizeMethod value=\"Cole2000\"/>\n\
    <majorMergerMassRatio value=\"0.25\"/>\n\
    <mergerRemnantSizeOrbitalEnergy value=\"1\"/>\n\
    \n\
    <!-- Spheroid options -->\n\
    <spheroidEnergeticOutflowMassRate value=\"1.0e-2\"/>\n\
    \n\
    <!-- Numerical tolerances -->\n\
    <odeToleranceAbsolute value=\"0.01\"/>\n\
    <odeToleranceRelative value=\"0.01\"/>\n\
    <diskMassToleranceAbsolute value=\"1.0e-6\"/>\n\
    <spheroidMassToleranceAbsolute value=\"1.0e-6\"/>\n\
    <timestepHostAbsolute value=\"1.0\"/>\n\
    <timestepHostRelative value=\"0.1\"/>\n\
    <timestepSimpleAbsolute value=\"1.0\"/>\n\
    <timestepSimpleRelative value=\"0.1\"/>\n\
    \n\
    <!-- Output options -->\n\
    <mergerTreeOutputReferences value=\"false\"/>\n\
    \n\
</parameters>\n\
";
    fs.close();
}

}
