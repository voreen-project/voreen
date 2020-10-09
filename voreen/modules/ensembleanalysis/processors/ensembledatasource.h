/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_ENSEMBLEDATASOURCE_H
#define VRN_ENSEMBLEDATASOURCE_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/ports/genericport.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/string/stringtableproperty.h"

#include "modules/plotting/properties/colormapproperty.h"

#include "../ports/ensembledatasetport.h"

namespace voreen {

/**
 * Loads an Ensemble of Volumes, organized into multiple members.
 */
class VRN_CORE_API EnsembleDataSource : public Processor {

    static const std::string SCALAR_FIELD_NAME; // Deprecated!
    static const std::string NAME_FIELD_NAME;
    static const std::string SIMULATED_TIME_NAME;
    static const std::string RUN_NAME;
    static const std::string FALLBACK_FIELD_NAME;

public:
    EnsembleDataSource();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "EnsembleDataSource";    }
    virtual std::string getCategory() const   { return "Input";                 }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }
    virtual bool usesExpensiveComputation() const { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Loads an Ensemble of Volumes, organized into multiple members.");
        ensemblePath_.setDescription("Expects a folder containing a separate folder for each ensemble member "
                                     "(aka member). Each folder must contain a single file (of any supported format) "
                                     "for each time step of the respective member. Each of those files may contain multiple "
                                     "volumes which will be interpreted as time steps, in lexicographic order. "
                                     "The volume files might contain time step meta information which do not need to "
                                     "match the files' names. ");
        loadingStrategy_.setDescription("Loading strategy has three options:<br>"
                                        "<strong>Manual</strong>: The ensemble is only loaded when pressing the load button.<br>"
                                        "<strong>Full</strong>: The entire ensemble is loaded fully from disk when the workspace is loaded<br>"
                                        "<strong>Lazy</strong>: The entire ensemble needs to be loaded once an all required meta data will be "
                                        "safed to and loaded from disk the next time the workspace is loaded");
    }

    void process();
    virtual void initialize();
    virtual void deinitialize();

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    void clearEnsembleDataset();
    void buildEnsembleDataset();
    void printEnsembleDataset();

    std::vector<std::unique_ptr<const VolumeBase>> volumes_;
    std::unique_ptr<EnsembleDataset> output_;

    FileDialogProperty ensemblePath_;
    StringOptionProperty loadingStrategy_;
    ButtonProperty loadDatasetButton_;
    ProgressProperty memberProgress_;
    ProgressProperty timeStepProgress_;
    StringTableProperty loadedMembers_;
    FileDialogProperty printEnsemble_;
    ColorMapProperty colorMap_;
    BoolProperty overrideTime_;
    StringProperty hash_;

    /// The structure of the ensemble data.
    EnsembleDatasetPort outport_;

    static const std::string loggerCat_;
};

} // namespace

#endif
