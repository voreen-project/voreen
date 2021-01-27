/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_LAMBDACALLBACK_H
#define VRN_LAMBDACALLBACK_H

#include "voreen/core/datastructures/callback/callback.h"
#include <functional>
namespace voreen {

#define ON_CHANGE_LAMBDA(PROPERTY,FUNCTION) PROPERTY.onChange(LambdaFunctionCallback(FUNCTION));
/**
 * This Callback can call a c++11 lambda function
 */
class LambdaFunctionCallback : public Callback {
public:

    LambdaFunctionCallback(std::function<void()> func)
        : func_(func)
    {
    }

    virtual ~LambdaFunctionCallback() {}
    virtual LambdaFunctionCallback* clone() const override { return new LambdaFunctionCallback(func_); }

    virtual void exec() {
        func_();
    }

private:
    std::function<void()> func_;
};

}   // namespace

#endif
