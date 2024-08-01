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

#ifndef VRN_CMHALODATASOURCE_H
#define VRN_CMHALODATASOURCE_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "../ports/cmhaloport.h"
#include "../datastructures/cmmergertree.h"
#include "boost/iostreams/device/mapped_file.hpp"
#include <vector>



namespace voreen{

/**
 * An input processor for halo/merger tree data. The data is expected to be in
 * a textual table format with one file per time step and one line per halo. See
 * the provided data of the 2015 SciVis contest.
 * On after reading in the data from the text files a binary cache file will be
 * created in voreens cache directory, identified by (a hash of) the file names
 * of the time step files.
 * Additional members for each halo can be added by exporting the current data
 * to galacticus, run galacticus using the generated parameter file and
 * reimporting the data. The additional data will be added to the cache file
 * and will be present in subsequent runs as long as the cache file is
 * present.
 */
class CMHaloDataSource : public Processor{
public:
    CMHaloDataSource();
    virtual ~CMHaloDataSource();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "CMHaloDataSource"; }
    virtual std::string getCategory() const     { return "Viscontest2019";   }
    virtual void setDescriptions()              {setDescription("Cosmology Halo Data Source");}
    virtual CodeState   getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }



protected:
    virtual void initialize();
    virtual void deinitialize();

    virtual void process();

private:

    /**
     * Close the currently open cache file.
     */
    void closeFile();

    /**
     * Save the tree currently held in memory (haloLists_) to a cache file in Voreen's cache folder
     */
    void saveTreeToFile();

    /**
     * Export the tree given as the haloLists parameter to the file given by exportToFilePath.
     * @param exportFilePath Export file path for galacticus input file (HDF5 format)
     * @param haloLists std::vector, that holds the halos to be exported
     */
    void exportToHDF5(const std::string& exportFilePath, std::vector<std::vector<CMHalo>>& haloLists) const;

    /**
     * Import the output data of galacticus and update the halos given as parameter hLists
     * @param importFilePath Galacticus output file (HDF5 format)
     * @param hLists std::vector of halos to be updated
     */
    bool importFromHDF5(const std::string& importFilePath, std::vector<std::vector<CMHalo>>& hlists) const;

    /**
     * Read a single time step of halos from the given file and write the halos into the vector given as a parameter.
     * @param fileName halo time step file
     * @param halos Output std::vector, the halos read in will be appended to
     */
    void readHaloListFile(const std::string& fileName, std::vector<CMHalo>& halos) const;

    /**
     * Read the halo lists from the given directory, create a cache file if necessary, memory map a new or existing
     * cache file and point mergerTree_ to it.
     */
    void readmergerTreeData();

    /**
     * Export the loaded merger tree for galacticus and create a parameter file with the name given by galacticusOutputFile_.
     */
    void exportToGalacticus();

    /**
     * Load the information provided by galacticus (galacticusOutputFile_) and update the merger tree.
     */
    void importFromGalacticus();

    /**
     * Will be called if the property that indicates the folder which holds the halo list files changes.
     * It updates the file name hash (haloFilesIdentifier_) accordingly.
     */
    void fileFolderChanged();

    /**
     * Mark the currently cached merger tree invalid.
     */
    void invalidateMergerTree();

    /**
     * Loads the current in memory merger tree if necessary and returns a reference to it.
     */
    std::vector<std::vector<CMHalo>>& getHaloLists();

    /**
     * Get the path to the current cache file.
     */
    std::string getCacheFilePath() const;

    /**
     * Gets the set of red shifts of the currently loaded merger tree as a space separated list.
     */
    std::string getRedshiftsString();

    /**
     * Write a galacticus configuration file (configFileName) which instructs galacticus to calculate additional
     * information on the exported merger tree (hdfFileName)
     */
    void writeGalacticusConfig(const std::string configFileName, const std::string& hdfFileName);

    /// Halo data output
    CMHaloPort outport_;
    /// Folder that contains the halo list files
    FileDialogProperty fileFolder_;
    /// Export destination for a HDF5 file to be used by galacticus
    FileDialogProperty galacticusInputFile_;
    /// Export file of galacticus that contains additional data for all halos in the merger tree
    FileDialogProperty galacticusOutputFile_;
    /// Button create a galacticus parameter and import file
    ButtonProperty galacticusExportButton_;
    /// Button to import data from galacticus export file
    ButtonProperty galacticusImportButton_;
    /// Memory mapped merger tree cache file
    boost::iostreams::mapped_file_source* mergerTreeFile_;
    /// Pointer to the current merger tree
    const CMMergerTree* mergerTree_;

    /// Indicates whether the halo files for the current fileFolder_ have been loaded, yet.
    bool haloFilesLoaded_;
    /// Indicates whether the current halo lists contain additional data calculated by galacticus
    bool haloListsContainGalacticusData_;
    /// hash of all halo list files to identify a merger tree. It is used to generate the cache file name.
    std::string haloFilesIdentifier_;
    /// In memory representation of the merger tree. Created on demand from the halo list files.
    std::vector<std::vector<CMHalo>> haloLists_;

    /// Indicates whether the current mergerTree_ is valid or has to be recalculated (and the cache file has to be overwritten).
    bool treeValid_;
};
}

#endif /* end of include guard: VRN_CMHALODATASOURCE_H */
