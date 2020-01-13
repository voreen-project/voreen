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

#ifndef VRN_TRANSFUNC1DGAUSSIAN_H
#define VRN_TRANSFUNC1DGAUSSIAN_H


#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"
#include "utils/transfuncmappingcurve.h"

#include "voreen/core/datastructures/volume/histogram.h"
#include "tgt/vector.h"
#include <queue>
#include "voreen/core/utils/glsl.h"

namespace voreen {

/**
 * One dimensional, piece-wise linear transfer function based on Gauss curves.
 *
 * Internally, it is represented by a one-dimensional RGBA texture of type GL_UNSIGNED_BYTE.
 */
class VRN_CORE_API TransFunc1DGaussian : public TransFunc1D {
public:
    /**
     * Constructor
     *
     * @param width desired width of the transfer function
     */
    TransFunc1DGaussian(int width = 1024);

    /**
     * Copy constructor.
     *
     * @param tf reference to the TransFuncIntensity instance, which is copied
     */
    TransFunc1DGaussian(const TransFunc1DGaussian& tf);

    /**
     * Destructor - deletes the curves of the transfer function
     */
    virtual ~TransFunc1DGaussian();

    virtual std::string getClassName() const { return "TransFunc1DGaussian";     }
    virtual TransFunc1DGaussian* create() const  { return new TransFunc1DGaussian(); }
    virtual TransFunc1DGaussian* clone() const;
    virtual void setMemberValuesFrom(const TransFuncBase* transfunc);
    virtual bool compareTo(const TransFuncBase& tf) const;
    virtual void reset();

    //--------------------------------------
    //  handle texture
    //--------------------------------------
protected:
    /** @override */
    virtual tgt::Vector4<GLfloat> getMappingForValueFloat(float x);
    /** @override */
    virtual tgt::Vector4<GLubyte> getMappingForValueUByte(float x);

    //--------------------------------------
    //  handle curves
    //--------------------------------------
public:
    /**
     * Returns the number of curves in this transfer function.
     *
     * @return the number of curves
     */
    int getNumCurves() const;

    /**
     * Returns the i-th curve. Curves are sorted by their activation state.
     * If a value outside of [0, getNumCurves()] is passed, it will be clamped to the appropriate values.
     *
     * @param i the i-th curve will be returned
     * @return the pointer to the appropriate curve
     */
    const TransFuncMappingCurve* getCurve(int i) const;

    TransFuncMappingCurve* getCurve(int i);

    /**
     * Returns all curves. Curves are sorted by their activation state.
     *
     * @return a vector containing all the curves
     */
    const std::vector<TransFuncMappingCurve*> getCurves() const;

    /**
    * Replaces the current curves by the passed ones.
    *
    * @param curves the new curves
    */
    void setCurves(std::vector<TransFuncMappingCurve*> curves);

    /**
    * Sorts the curves in this gaussian transfer function in the following order:
    * deactivated curves, active curves, selected curve
    *
    * @param selectedCurve the currently selected curve
    */
    void sortCurves(TransFuncMappingCurve* selectedCurve);

    /**
     * Adds a curve to the property mapping function.
     * It will be automatically inserted into the correct position.
     *
     * @param curve the curve to be added
     */
    void addCurve(TransFuncMappingCurve* curve);

    /**
     * Updates a curve within the property mapping function.
     *
     * @param curve the curve to be updated
     */
    void updateCurve(TransFuncMappingCurve* curve);

    /**
     * Removes the given curve from the transfer function.
     * Also deletes the passed curve.
     *
     * @param curve the curve that will be removed.
     */
    void removeCurve(TransFuncMappingCurve* curve);

    /**
     * Removes all curves from the transfer function.
     */
    void clearCurves();

    /**
     * This method returns whether the mapping function is empty.
     * i.e., it contains no mapping curves.
     *
     * @return Is the mapping function empty?
     */
    bool isEmpty() const;
protected:
    /**
     * Calculates the average for the segment [segStart,segEnd).
     *
     * @param segStart the start of the segment
     * @param segEnd the end of the segment; the value itself is not part of the calculation
     * @return a vector containing the average value for the segment
     */
    tgt::vec4 getMeanValue(float segStart, float segEnd) const;

    //--------------------------------------
    //  function helper
    //--------------------------------------
public:
    /**
     * Creates a default function.
     * Generates a single white curve at intensity 0.5.
     */
    void setToStandardFunc();

    /*
     * Calculates a gaussian transfer function with two curves based on the
     * automatic treshold algorithm by Otsu
     */
    void autoTreshold(Histogram1D* histogram);

    //--------------------------------------
    //  load and save
    //--------------------------------------
public:
    virtual const std::vector<std::string> getLoadFileFormats() const;
    virtual const std::vector<std::string> getSaveFileFormats() const;

        /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);

    /**
     * The central entry point for loading a transfer function. The file extension must be tfg.
     *
     * @param filename the name of the file, which should be opened
     * @return true if loading succeeds, false otherwise
     */
    virtual bool load(const std::string& filename);

    /**
     * Saves the transfer function to a tfg-file. Any data in the file will be overwritten.
     *
     * @param filename the name of the file the transfer function will be saved to
     * @return true, if the operation was successfull, false otherwise
     */
    virtual bool save(const std::string& filename) const;
protected:
    /**
     * Saves transfer function to a XML file. The extension of the file is tfi.
     *
     * @param filename the name of the file the transfer function will be saved to
     * @return true, if the operation was successful, false otherwise
     */
    bool saveTfg(const std::string& filename) const;

        /**
     * Loads a transfer function from an file with ending tfi.
     *
     * @param filename the name of the file, which should be opened
     * @return true if loading succeeds, false otherwise
     */
    bool loadTfg(const std::string& filename);


protected:
    std::vector<TransFuncMappingCurve*> curves_; ///< internal representation of the transfer function as a set of curves

    static const std::string loggerCat_; ///< logger category
};

///------------------------------------------------------------------------------------------------

/**
 * Metadata encapsulating TransFunc1DGaussian.
 */
class VRN_CORE_API TransFunc1DGaussianMetaData : public TransFuncMetaDataGeneric<TransFunc1DGaussian> {
public:
    TransFunc1DGaussianMetaData();
    TransFunc1DGaussianMetaData(TransFunc1DGaussian* transfunc);
    TransFunc1DGaussianMetaData(const std::vector<TransFunc1DGaussian*>& transfunc);
    virtual MetaDataBase* clone() const;
    virtual MetaDataBase* create() const     { return new TransFunc1DGaussianMetaData(); }
    virtual std::string getClassName() const { return "TransFunc1DGaussianMetaData";     }
    virtual std::string toString() const {return "TransFunc1DGaussianMetaData";}
};


} // namespace voreen

#endif // VRN_TRANSFUNC1DGAUSSIAN_H
