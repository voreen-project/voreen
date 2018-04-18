/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#ifndef VRN_STREAMLINESOURCE_H
#define VRN_STREAMLINESOURCE_H

#include "voreen/core/processors/processor.h"
#include "../../ports/streamlinelistport.h"
#include "voreen/core/ports/volumeport.h"

#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/string/stringtableproperty.h"

namespace voreen {

class StreamlineList;

/**
 * Used to load streamlines previously stored with the StreamlineSave processor.
 * The load routine is implemented in the Streamline class itself.
 *
 * @see StreamlineSave
 * @see Streamline
 */
class VRN_CORE_API StreamlineSource : public Processor {

public:
    StreamlineSource();

    virtual Processor* create() const override          { return new StreamlineSource(); }
    virtual std::string getClassName() const override   { return "StreamlineSource";     }
    virtual std::string getCategory() const override    { return "Input";                }
    virtual CodeState getCodeState() const override     { return CODE_STATE_TESTING;     }

protected:
    virtual void setDescriptions() override {
        setDescription("Loads streamlines previously stored with the <i>StreamlineSave</i> processor (VSD-Format). The loaded streamlines can be visualized or modified " \
                        "by other processors of the <i>Flowreen</i> module.");
        loadOptionProp_.setDescription("Switch between loading a single StreamlineList and an entire directory.");
        //single
        singleFileProp_.setDescription("Single VSD file, which should be loaded.");
        //directory
        directoryFileProp_.setDescription("All VSD files of the selected directory are loaded.");
        tableProp_.setDescription("This property shows all loaded StreamlineLists. The currently selected list is highlighted.");
        currentListProp_.setDescription("The currently selected StreamlineList.");
        //ports
        streamlineOutport_.setDescription("Port containing the currently selected/loaded Streamlinelist.");
        volumeOutport_.setDescription("Port containing the magnitude volume or NULL, if no volume is contained in the vsd file.");
    }
    virtual void process() override;
    /** Triggers a auto-save if the file property has changed values. */
    virtual void invalidate(int inv) override;
    /** Volume outport connection should be optional.*/
    virtual bool isReady() const override;


    /** Determines if a single file or a directory should be selected. */
    enum LoadOption {
        LOAD_SINGLE,
        LOAD_DIRECTORY
    };

    /** Function used to load a single streamline file. */
    void loadSingleStreamlines();
    /** Function used to load a single streamline file. */
    void loadDirectoryStreamlines();

    //--------------
    //  Helper
    //--------------
    /** Returns the loaded file or null.*/
    std::pair<StreamlineList*,Volume*> loadVSDFile(const std::string& absPath);

    //--------------
    //  Callbacks
    //--------------
    /** Triggered by loadOptionProp_. Changes property visibilities. */
    void loadOptionOnChange();
    /** Triggered by currentListProp_. Forces an output update. */
    void currentListOnChange();
    /** Triggered by tableProp_. Updates the output via currentListOnChange. */
    void tableSelectionOnChange();

    //--------------
    //  Member
    //--------------
private:
    //ports
    StreamlineListPort streamlineOutport_;      ///< port containing the loaded streamlines
    VolumePort volumeOutport_;                  ///< port containing the refenrece magnitude image, if present

    //properties
    OptionProperty<LoadOption> loadOptionProp_; ///< Property defines, if a single file or a directory should be loaded
        //single
    FileDialogProperty singleFileProp_;         ///< name of a single file to load
        //directory
    FileDialogProperty directoryFileProp_;      ///< name of the directory containing all files
    StringTableProperty tableProp_;             ///< table containing all loaded files
    IntProperty currentListProp_;               ///< currently selected list

    bool forceStreamlineLoadFlag_;              ///< if true, the streamlines will be loaded in process
    bool updateOutputFlag_;                     ///< updates the output on currentListProp_ changes
    std::vector<std::pair<StreamlineList*,Volume*> > currentlyLoadedListsAndVolumes_; ///< vector containing the currently loaded lists/volumes in directory mode

    bool firstDirectoryLoad_;                   ///< used to deserialie the last set list correclty

    static const std::string loggerCat_;        ///< Meow
};

} // namespace

#endif
