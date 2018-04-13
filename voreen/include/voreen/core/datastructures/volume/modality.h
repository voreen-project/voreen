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

#ifndef VRN_MODALITY_H
#define VRN_MODALITY_H

#include <string>
#include <iostream>
#include <vector>

namespace voreen {

class Modality {
public:
    Modality(const std::string& name = "");

    /// Two modalities are equal if their names are equal
    bool operator== (const Modality& m) const { return m.name_ == name_; }
    bool operator!= (const Modality& m) const { return m.name_ != name_; }
    friend std::ostream& operator<<(std::ostream& os, const Modality& m);

    std::string getName() const { return name_; }
    bool isUnknown() const { return name_ == MODALITY_UNKNOWN.name_; }

    static const std::vector<Modality*> getModalities() { return modalities(); }
    static const std::vector<std::string>& getModalityNames() { return modalityNames_(); }

    // all known modalities
    static const Modality MODALITY_UNKNOWN;
    static const Modality MODALITY_ANY;
    static const Modality MODALITY_CT;
    static const Modality MODALITY_PET;
    static const Modality MODALITY_MR;
    static const Modality MODALITY_US;
    static const Modality MODALITY_SEGMENTATION;
    static const Modality MODALITY_MASKING;
    static const Modality MODALITY_LENSEVOLUME;
    static const Modality MODALITY_AMBIENTOCCLUSION;
    static const Modality MODALITY_DYNAMICAMBIENTOCCLUSION;
    static const Modality MODALITY_DIRECTION_X;
    static const Modality MODALITY_DIRECTION_Y;
    static const Modality MODALITY_DIRECTION_Z;
    static const Modality MODALITY_DIRECTIONS;
    static const Modality MODALITY_NORMALS;
    static const Modality MODALITY_GRADIENTS;
    static const Modality MODALITY_GRADIENT_MAGNITUDES;
    static const Modality MODALITY_2ND_DERIVATIVES;
    static const Modality MODALITY_FLOW;
    static const Modality MODALITY_INDEX_VOLUME;
    static const Modality MODALITY_BRICKED_VOLUME;
    static const Modality MODALITY_EEP_VOLUME;

protected:
    std::string name_;

    // This indirection is done to prevent the "static initialization order fiasco"
    // see http://www.parashift.com/c++-faq-lite/ctors.html (10.13)
    static std::vector<std::string>& modalityNames_() {
        static std::vector<std::string> modalityNames;
        return modalityNames;
    }

    static std::vector<Modality*>& modalities() {
        static std::vector<Modality*> modalities;
        return modalities;
    }
};

inline std::ostream& operator<<(std::ostream& os, const Modality& m) {
    return os << m.name_;
}

} // namespace voreen

#endif // VRN_MODALITY_H
