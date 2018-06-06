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

#ifndef VRN_PPSYMBOL_H
#define VRN_PPSYMBOL_H

#include "voreen/core/utils/GLSLparser/symbol.h"

#include "voreen/core/utils/GLSLparser/preprocessor/ppparsetreenode.h"
#include <vector>

namespace voreen {

namespace glslparser {

class PreprocessorSymbol : public Symbol {
public:
    PreprocessorSymbol(const std::string& id, const bool isFunction)
        : Symbol(id),
        defined_(true),
        isFunction_(isFunction),
        body_(0)
    {
    }

    PreprocessorSymbol(const std::string& id, const bool isFunction, TokenList* const body)
        : Symbol(id),
        defined_(true),
        isFunction_(isFunction),
        body_(body)
    {
    }

    PreprocessorSymbol(const std::string& id, const bool isFunction, TokenList* const body,
                       const std::vector<std::string>& formals)
        : Symbol(id),
        defined_(true),
        isFunction_(isFunction),
        body_(body),
        formals_(formals)
    {
    }

    ~PreprocessorSymbol() {
        delete body_;
    }

    bool isDefined() const { return defined_; }
    bool isFunction() const { return isFunction_; }
    TokenList* getBody() const { return body_; }
    const std::vector<std::string>& getFormals() const { return formals_; }
    size_t getNumFormals() const { return formals_.size(); }

    void undefine() { defined_ = false; }

protected:
    bool defined_;
    bool isFunction_;
    TokenList* const body_;
    std::vector<std::string> formals_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
