/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_LINKEVALUATORTRANSFERFUNCTIONPROPERTIES_H
#define VRN_LINKEVALUATORTRANSFERFUNCTIONPROPERTIES_H

#include "voreen/core/properties/link/linkevaluatorbase.h"
#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"

namespace voreen {

class LinkEvaluatorTransFunc1DPropertyBase : public LinkEvaluatorBase {
    virtual void eval(Property* src, Property* dst);
    virtual bool arePropertiesLinkable(const Property* p1, const Property* p2) const;

    virtual float getTFProperty(TransFunc1D& tf) = 0;
    virtual void setTFProperty(TransFunc1D& tf, float value) = 0;
};

class LinkEvaluatorTransFunc1DThresholdMin : public LinkEvaluatorTransFunc1DPropertyBase {
public:
    virtual float getTFProperty(TransFunc1D& tf);
    virtual void setTFProperty(TransFunc1D& tf, float value);

    virtual std::string getClassName()  const { return "LinkEvaluatorTransFunc1DThresholdMin"; }
    virtual std::string getGuiName()    const { return "Transfunc1D (minimum threshold) <-> Float"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorTransFunc1DThresholdMin(); }
};

class LinkEvaluatorTransFunc1DThresholdMax : public LinkEvaluatorTransFunc1DPropertyBase {
public:
    virtual float getTFProperty(TransFunc1D& tf);
    virtual void setTFProperty(TransFunc1D& tf, float value);

    virtual std::string getClassName()  const { return "LinkEvaluatorTransFunc1DThresholdMax"; }
    virtual std::string getGuiName()    const { return "Transfunc1D (maximum threshold) <-> Float"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorTransFunc1DThresholdMax(); }
};

class LinkEvaluatorTransFunc1DDomainMin : public LinkEvaluatorTransFunc1DPropertyBase {
public:
    virtual float getTFProperty(TransFunc1D& tf);
    virtual void setTFProperty(TransFunc1D& tf, float value);

    virtual std::string getClassName()  const { return "LinkEvaluatorTransFunc1DDomainMin"; }
    virtual std::string getGuiName()    const { return "Transfunc1D (minimum domain) <-> Float"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorTransFunc1DDomainMin(); }
};

class LinkEvaluatorTransFunc1DDomainMax : public LinkEvaluatorTransFunc1DPropertyBase {
public:
    virtual float getTFProperty(TransFunc1D& tf);
    virtual void setTFProperty(TransFunc1D& tf, float value);

    virtual std::string getClassName()  const { return "LinkEvaluatorTransFunc1DDomainMax"; }
    virtual std::string getGuiName()    const { return "Transfunc1D (maximum domain) <-> Float"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorTransFunc1DDomainMax(); }
};

class LinkEvaluatorTransFunc1DGamma : public LinkEvaluatorTransFunc1DPropertyBase {
public:
    virtual float getTFProperty(TransFunc1D& tf);
    virtual void setTFProperty(TransFunc1D& tf, float value);

    virtual std::string getClassName()  const { return "LinkEvaluatorTransFunc1DGamma"; }
    virtual std::string getGuiName()    const { return "Transfunc1D (gamma) <-> Float"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorTransFunc1DGamma(); }
};

}
#endif
