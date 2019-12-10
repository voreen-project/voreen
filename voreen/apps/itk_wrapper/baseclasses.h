/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef BASECLASSES_H
#define BASECLASSES_H


namespace voreen {

/***********************************************************************************
 * Class Argument
 * An Argument is part of a (ITK-)Filter which is represented as a property
 * in Voreen.
 ***********************************************************************************/
class Argument : public Serializable {
public:
    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);

    /**
     * Returns the voreen-propertyname of the argument.
     * (e.g. argument.name_ = "ThresholdValue", propertyname = "thresholdValue_")
     *
     * @return propertyname, name of the argument's voreen-property
    */
    std::string getPropertyname();

    /**
     * Returns the property-initialization of the argument for the filter-.cpp-file.
     *
     * @return result, property-initialization of the argument
    */
    std::string getProperty();

    /**
     * Returns the ITK-filter->SetArgument part for an input-argument
     * or the ITK-filter->GetArgument part for an output_argument of the filter.
     *
     * @return result, filter->SetArgument part or filter->GetArgument part of the argument
    */
    std::string getFilterSetOrGet();

    /**
     * Returns the name of the argument
     *
     * @return name_, name of the argument
    */
    std::string getName();

    /**
     * Returns the type of the argument.
     *
     * @return type_, type of the argument
    */
    std::string getType();

    /**
     * Returns whether the argument is an input- or an output-argument.
     *
     * @return inputArgument_, true if the argument is an input-argument and false if it is an output-argument
    */
    bool isInputArgument();

    /**
     * Returns the default value of the argument.
     *
     * @return defaultVaue_, default value of the argument
    */
    std::string getDefaultValue();

    /**
     * Returns the minimum value of the argument
     *
     * @return minVaue_, minimum value of the argument
    */
    std::string getMinValue();

    /**
     * Returns the maximum value of the argument.
     *
     * @return maxVaue_, maximum value of the argument
    */
    std::string getMaxValue();

private:
    /**
     * Returns the property-ID of the argument.
     * (This is the first argument of each property.)
     *
     * @return propertyID, property-ID of  of the argument
    */
    std::string getPropertyID();

    /**
     * Returns the a temporary name which is needed for some special argument-types.
     * (e.g. ArrayType)
     *
     * @return tempName, temporary name for some special argument-types
    */
    std::string getTempName();

    std::string name_;          ///< Name of the argument
    std::string type_;          ///< Type of the argument
    bool inputArgument_;        ///< Boolean defining whether the argument is an input- or an output-argument
    std::string defaultValue_;  ///< Default value of the argument
    std::string minValue_;      ///< Minimum value of the argument
    std::string maxValue_;      ///< Maximum value of the argument

};

/***********************************************************************************
 * Class Ports
 * A Port is used to define the input- or output-port of a (ITK-)Filter in Voreen.
 ***********************************************************************************/
class Ports : public Serializable {
public:
    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);

    /**
     * Returns an ITK-voxel-type for a Voreen-volume-type.
     * (e.g. uint8_t for volume-type UInt8)
     *
     * @param index, position of the volume-type in possibleTypes_
     * @return typeName, voxel-type for the volume-type at the index-position
     */
    std::string getITKTypename(int index);

    /**
     * Returns the name of the port.
     *
     * @return name_, name of the port
    */
    std::string getName();

    /**
     * Sets the name of the port.
     *
     * @param name
    */
    void setName(std::string name);

    /**
     * Returns whether the port is a template-argument of the filter or not.
     *
     * @return filterTemplate_, true if the port is a template-argument of the filter and false if not
    */
    bool isFilterTemplate();

    /**
     * Sets filterTemplate_ to bool filterTemplate
     *
     * @param filterTemplate
    */
    void setFilterTemplate(bool filterTemplate);

    /**
     * Returns whether the the name of the port is used as filter-input-setter or not.
     * (Usually the input of a filter is set by filter->SetInput(...), but for some special
     * inputs the name of the inport is used as setter, e.g. filter->SetMaskImage(...))
     *
     * @return nameIsSetter_, true if the name of the inport is the filter-input-setter and false if not
    */
    bool isNameIsSetter();

    /**
     * Sets nameIsSetter_ to bool nameIsSetter.
     *
     * @param nameIsSetter
    */
    void setNameIsSetter(bool nameIsSetter);

    /// Should the output on this port get the RealWorldMapping from the inport?
    bool transferRWM();
    void setTransferRMW(bool tRWM);

    /// Should the output on this port get the transformation matrix from the inport?
    bool transferTransformation();
    void setTransferTransformation(bool tTrafo);

    /**
     * Returns the possible volume-types for a port.
     *
     * @return possibleTypes_, possible volume-types for a port
    */
    std::vector<std::string> getPossibleTypes();

    /**
     * Sets the possibleTypes_ for a port to types.
     *
     * @param types, possible volume-types for a port
    */
    void setPossibleTypes(std::vector<std::string> types);

    /**
     * Adds a type to possibleTypes_ of a port.
     *
     * @param type, possible volume-type to add to possibleTypes_ of a port
    */
    void addPossibleType(std::string type);

    /**
     * Deletes a metaType (e.g. "ScalarTypes") from possibleTypes_ of a port.
     * (This is done after adding all types of a metaType to possibleTypes_ of a port.)
     *
     * @param position, position of the metaType in possibleTypes_
    */
    void deleteMetaType(size_t position);

private:
    std::string name_;                          ///< Name of a port
    bool filterTemplate_;                       ///< Boolean defining whether the port is a template-argument of the filter or not
    bool nameIsSetter_;                         ///< Boolean defining whether the name of the port is the filter-input-setter of or not
    bool transferRWM_;                          ///< Boolean defining whether the output on this port should get the RealWorldMapping from the inport
    bool transferTransformation_;               ///< Boolean defining whether the output on this port should get the physicalToWorldMatrix from the inport
    std::vector<std::string> possibleTypes_;    ///< Vector containing possible volume-types for the port
};

/***********************************************************************************
 * Class Filter
 * A Filter provides methods to create a Voreen-processor (.cpp- and .h-file)
 * which represents an ITK-Filter.
 ***********************************************************************************/
class Filter : public Serializable {
public:
    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);

    /**
     * Initializes the ports of a filter.
     * This contains:
     * - If a filter has no port specified in the xml-file, it becomes a standard filter
     *   with one inport and one outport with scalar-volume-types as possibleTypes_.
     * - If the possibleTypes_-varible of a port of the filter contains a metaType, it will
     *   be replaced by all types belonging to the metaType.
     * - The variables scalarVolumeTypesPossible_, vec2VolumeTypesPossible_, vec3VolumeTypesPossible_
     *   and vec4VolumeTypesPossible_ will be set after checking all possible types of each port.
     *
     * @param integerScalarTypes, vector containing all types of the metaType "IntegerScalar"
     * @param realScalarTypes, vector containing all types of the metaType "RealScalar"
     * @param scalarTypes, vector containing all types of the metaType "Scalar"
     * @param integerVectorTypes, vector containing all types of the metaType "IntegerVector"
     * @param realVectorTypes, vector containing all types of the metaType "RealVector"
     * @param vectorTypes, vector containing all types of the metaType "Vector"
     */
    void initializeFilterPorts(std::vector<std::string> integerScalarTypes
        , std::vector<std::string> realScalarTypes
        , std::vector<std::string> scalarTypes
        , std::vector<std::string> integerVectorTypes
        , std::vector<std::string> realVectorTypes
        , std::vector<std::string> vectorTypes);

    /**
     * Creates the .cpp-file of a filter
     *
     * @param moduleName, name of the filters module
     */
    void createCPP(std::string moduleName);

    /**
     * Creates the .h-file of a filter
     *
     * @param moduleName, name of the filters module
     * @param group, group of the filters module
     */
    void createHeaderfile(std::string moduleName, std::string group);

    /**
     * Returns the filename of a filter (all letters to lower-case).
     *
     * @return filename
     */
    std::string getFilename();

    /**
     * Returns the header-guard of a filter (all letters to upper-case).
     *
     * @return headerGuard
     */
    std::string getHeaderGuard();

    /**
     * Returns the path of a special template-file.
     *
     * @param templateFilename
     * @return path, path of the template with the filename templateFilename
     */
    std::string getTemplatePath(std::string templateFilename);

    /**
     * Returns the path of the filter in an ITK-module-directory.
     *
     * @param moduleName
     * @param filterFilename
     * @return path, path of the filter in "moduleName"-directory
     */
    std::string getProcessorOutputPath(std::string moduleName, std::string filterFilename);

    /**
     * Returns the name of the filter.
     *
     * @return name_, name of the port
    */
    std::string getName();

    /**
     * Returns whether the filter is enabled or not.
     *
     * @return enabled_, true if the filter is enabled and false if it is disabled
    */
    bool isEnabled();

    /**
     * Returns the disable-description of the filter.
     *
     * @return description_, description explaining why the filter is disabled
    */
    std::string getDescription();

    /**
     * Returns whether the filter can be auto-generated by the wrapper or not.
     *
     * @return autoGenerated_, true if the filter can be auto-generated by the wrapper else false
    */
    bool isAutoGenerated();

    /**
     * Returns the code-state of the filter.
     *
     * @return codeState_, code-state of the filter
    */
    std::string getCodeState();

    /**
     * Returns the inports-vector of the filter.
     *
     * @return inports_, vector containing the inports of the filter
    */
    std::vector<Ports> getInports();

    /**
     * Returns the outports-vector of the filter.
     *
     * @return outports_, vector containing the outports of the filter
    */
    std::vector<Ports> getOutports();

    /**
     * Returns the arguments-vector of the filter.
     *
     * @return arguments_, vector containing the arguments of the filter
    */
    std::vector<Argument> getArguments();

private:
    /**
     * Creates a filter-method of the type "methodType".
     * The filter-method contains all filter-settings and if an outport exists the
     * outport-setting.
     *
     * @param methodType, type of the filter-method (e.g. "Vec2", empty for standard-filter-method)
     * @param kernel, true if the filter has a kernel as argument else false
     * @return filterMethod, string containing the filter-method
     */
    std::string getFilterMethod(std::string methodType, bool kernel);

    /**
     * Creates the processing-methods of a filter.
     * These proove the inport-volume-types and call the filter-methods to process.
     *
     * @return processingMethods, string containing the processing-methods
     */
    std::string getProcessingMethods();

    /**
     * Returns the name of the process-function.
     *
     * @return processFunctionName
     */
    std::string getProcessFunctionName();

    std::string name_;                  ///< Name of a filter
    bool enabled_;                      ///< Boolean whether the filter is enabled or disabled
    std::string description_;           ///< Has two meanings:
                                        ///< If a filter is disabled the description describes why it is disabled.
                                        ///< If the filter is enabled the description can be used to describe how the filter works.
    bool autoGenerated_;                ///< Boolean whether the filter can be auto-generated by the wrapper or not
    std::string codeState_;             ///< Code-state of a filter
    std::vector<Ports> inports_;        ///< Vector containing the inports of a filter
    std::vector<Ports> outports_;       ///< Vector containing the outports of a filter
    std::vector<Argument> arguments_;   ///< Vector containing the arguments of a filter

    bool scalarVolumeTypesPossible_;    ///< Boolean whether scalar-volume-types are possible or not
    bool vec2VolumeTypesPossible_;      ///< Boolean whether vec2-volume-types are possible or not
    bool vec3VolumeTypesPossible_;      ///< Boolean whether vec3-volume-types are possible or not
    bool vec4VolumeTypesPossible_;      ///< Boolean whether vec4-volume-types are possible or not
};

/***********************************************************************************
 * Class ITK_Module
 * An ITK_Module represents a module of ITK and contains all filters belonging
 * to this module.
 ***********************************************************************************/
class ITK_Module : public Serializable {
public:
    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);

    /**
     * Returns the name of the ITK-module.
     *
     * @return name_, name of the ITK-module
    */
    std::string getName();

    /**
     * Returns the group of the ITK-module.
     *
     * @return group_, group of the ITK-module
    */
    std::string getGroup();

    /**
     * Returns the filters-vector of the ITK-module.
     *
     * @return filters_, vector containing the filters of the ITK-module
    */
    std::vector<Filter> getFilters();

private:
    std::string name_;              ///< Name of a module
    std::string group_;             ///< Group of a module (e.g. "Filtering")
    std::vector<Filter> filters_;   ///< Vector containing the filters of a module
};

} //namespace voreen

#endif // BASECLASSES_H
