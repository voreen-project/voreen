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

#ifndef VRN_COMMANDLINEPARSER_H
#define VRN_COMMANDLINEPARSER_H

#include "voreen/core/voreencoreapi.h"
#include "voreen/core/utils/exception.h"
#include "tgt/assert.h"

#ifdef WIN32
#pragma warning(disable:4267) //< conversion from 'size_t' to 'unsigned int'
#endif
#include <boost/program_options.hpp>
#ifdef WIN32
#pragma warning(default:4267)
#endif

namespace po = boost::program_options;

namespace voreen {

/**
 * Command line parser for extracting program arguments,
 * based on the boost::program_options library.
 *
 * The parser supports typed, single-valued and multi-valued program options,
 * which are grouped into three categories: required, main, additional
 * The extracted option values can be written to reference parameters
 * or queried after the command line has been processed.
 * Additionally, all command line options can be read from a config file.
 *
 * @see VoreenApplication
 */
class VRN_CORE_API CommandLineParser {
public:

    /// Type of a program argument.
    enum OptionType {
        RequiredOption = 1,     ///< Required argument, i.e., the parser triggers an exception, if it is missing
        MainOption = 2,         ///< Optional argument that is to be listed as main option in the usage message
        AdditionalOption = 4,   ///< Optional argument that is to be listed as additional option in the usage message
        AllTypes = RequiredOption | MainOption | AdditionalOption
    };

    /**
     * Constructor.
     *
     * @param executableName The name of the program executable, used in the usage message
     * @param description Description of the program, added to the usage message
     * @param cmdStyle style of the command line
     * @param configFile path to the config file to read in addition to the command line (INI file syntax).
     *  If an empty string is passed, the config file is disabled. If the passed config file
     *  is not found, it is ignored. Command line options have preference over config file options.
     */
    CommandLineParser(const std::string& executableName, const std::string& description = "",
        po::command_line_style::style_t cmdStyle = po::command_line_style::default_style,
        const std::string& configFile = "");

    ~CommandLineParser();

    /**
     * Sets the style to use when parsing the command line.
     *
     * @note Needs to be called before execute() to have an effect.
     */
    void setCommandLineStyle(po::command_line_style::style_t cmdStyle);

    /**
     * Sets the path to the config file to read (INI file syntax) in addition to the command line.
     * If an empty string is passed, the config file is disabled. If the passed config file
     * is not found, it is ignored. Command line options have preference over config file options.
     *
     * @note Needs to be called before execute() to have an effect.
     */
    void setConfigFile(const std::string& configFile);

    /**
     * Sets the command line as passed to main().
     *
     * @param argc The number of arguments
     * @param argv The argument array
     */
    void setCommandLine(int argc, char** argv);

    /**
     * Parses command line and config file and evaluates the arguments.
     * Command line options have preference over config file options.
     *
     * @note Call this method after the command line has been assigned
     *  and the options have been added.
     *
     * @throw VoreenException if unknown or malformed arguments have been passed
     *  or if a required argument is missing
     */
    void execute();

    /**
     * Adds a typed, single-valued program option (--argument=value / --argument value).
     *
     * @tparam T Type of the program option. Note: the istream& operator>> has to be defined for it.
     *
     * @param name Name of the program option without dashes. Use "argument,a" for allowing short-form,
     *      e.g. "-a value" instead of "--argument value".
     * @param type Type type of the program option
     * @param description Description of the option to be shown in the usage message
     */
    template<class T>
    void addOption(const std::string& name, OptionType type = MainOption, const std::string& description = "");

    /**
     * Adds a typed, single-valued program option (--argument=value / --argument value).
     *
     * @tparam T Type of the program option. Note: the istream& operator>> has to be defined for it.
     *
     * @param name Name of the program option without dashes. Use "argument,a" for allowing short-form,
     *      e.g. "-a value" instead of "--argument value".
     * @param type Type type of the program option
     * @param description Description of the option to be shown in the usage message
     * @param defaultValue default value that is used, if the argument has not been specified by the user
     * @param defaultValueText textual representation of the default value
     */
    template<class T>
    void addOption(const std::string& name, OptionType type, const std::string& description,
        const T& defaultValue, const std::string& defaultValueText = "");

    /**
     * Adds a typed, single-valued program option (--argument=value / --argument value).
     *
     * @tparam T Type of the program option. Note: the istream& operator>> has to be defined for it.
     *
     * @param name Name of the program option without dashes. Use "argument,a" for allowing short-form,
     *      e.g. "-a value" instead of "--argument value".
     * @param out reference parameter the extracted argument value will be written to
     * @param type Type type of the program option
     * @param description Description of the option to be shown in the usage message
     */
    template<class T>
    void addOption(const std::string& name, T& out, OptionType type = MainOption, const std::string& description = "");

    /**
     * Adds a typed, single-valued program option (--argument=value / --argument value).
     *
     * @tparam T Type of the program option. Note: the istream& operator>> has to be defined for it.
     *
     * @param name Name of the program option without dashes. Use "argument,a" for allowing short-form,
     *      e.g. "-a value" instead of "--argument value".
     * @param out reference parameter the extracted argument value will be written to
     * @param type Type type of the program option
     * @param description Description of the option to be shown in the usage message
     * @param defaultValue default value that is used, if the argument has not been specified by the user
     * @param defaultValueText textual representation of the default value
     */
    template<class T>
    void addOption(const std::string& name, T& out, OptionType type, const std::string& description,
        const T& defaultValue, const std::string& defaultValueText = "");

    /**
     * Adds a boolean program option that does not expect a value (--argument).
     * If the argument is present in the command line, the option's value is true, otherwise false.
     *
     * @param name Name of the program option without dashes. Use "argument,a" for allowing short-form,
     *      e.g. "-a" instead of "--argument".
     * @param type Type type of the program option
     * @param description Description of the option to be shown in the usage message
     */
    void addFlagOption(const std::string& name, OptionType type = MainOption, const std::string& description = "");

    /**
     * Adds a boolean program option that does not expect a value (--argument).
     * If the argument is present in the command line, the option's value is true, otherwise false.
     *
     * @param name Name of the program option without dashes. Use "argument,a" for allowing short-form,
     *      e.g. "-a" instead of "--argument".
     * @param out reference parameter the flag value will be written to
     * @param type Type type of the program option
     * @param description Description of the option to be shown in the usage message
     */
    void addFlagOption(const std::string& name, bool& out, OptionType type = MainOption, const std::string& description = "");

    /**
     * Adds a typed multi-valued program option (--argument value1 value2 ...)
     *
     * @tparam T Type of the program option's values. Note: the istream& operator>> has to be defined for it
     *
     * @param name Name of the program option. For short-form, e.g. "-a" instead of "--argument", use the format: "argument,a"
     * @param type Type type of the program option
     * @param description Description of the option to be shown in the usage message
     */
    template<class T>
    void addMultiOption(const std::string& name, OptionType type = MainOption, const std::string& description = "");

    /**
     * Adds a typed multi-valued program option (--argument value1 value2 ...)
     *
     * @tparam T Type of the program option's values. Note: the istream& operator>> has to be defined for it
     *
     * @param name Name of the program option. For short-form, e.g. "-a" instead of "--argument", use the format: "argument,a"
     * @param out vector the extracted argument values will be written to
     * @param type Type type of the program option
     * @param description Description of the option to be shown in the usage message
     */
    template<class T>
    void addMultiOption(const std::string& name, std::vector<T>& out, OptionType type = MainOption, const std::string& description = "");

    /**
     * Retrieves the value of a single-valued or flag option extracted from the command line.
     *
     * @note Only valid after the command line has been processed (see execute()).
     */
    template<class T>
    void getOptionValue(const std::string& name, T& out) const;

    /**
     * Retrieves the values of a multi-valued option extracted from the command line.
     *
     * @note Only valid after the command line has been processed (see execute()).
     */
    template<class T>
    void getMultiOptionValues(const std::string& name, std::vector<T>& out) const;

    /**
     * Returns whether an option has been specified on the command line.
     *
     * @note Only valid after the command line has been processed (see execute()).
     */
    bool isOptionSet(const std::string& name) const;

    /**
     * Returns the first command line argument containing the path and the executable name.
     */
    std::string getProgramPath() const;

    /**
     * Returns all program options of the specified type that have been added to the parser,
     * regardless of whether they are present in the command line.
     */
    const po::options_description* getOptions(OptionType types = AllTypes) const;

    /**
     * Returns a usage message for the application consisting of the executable name, the program description
     * and the parameter description.
     *
     * @param optionTypes the types of program options to include in the message
     * @param programDescription if true, the description is included in the message
     */
    std::string getUsageString(OptionType optionTypes = AllTypes, bool programDescription = true) const;

protected:
    /// The name of the program executable, which is added to the usage string
    std::string executableName_;

    /// The description of the program, which is added to the usage string
    std::string programDescription_;

    /// Style to use when parsing the command line
    po::command_line_style::style_t cmdStyle_;

    /// Path to the config file to read.
    std::string configFile_;

    /// Program executable path as extracted from the command line
    std::string programPath_;

    int argc_;      ///< number of command line arguments
    char** argv_;   ///< raw command line

    po::options_description requiredOptions_;       ///< added required options
    po::options_description mainOptions_;           ///< added main options
    po::options_description additionalOptions_;     ///< added additional options

    /// temp variable used for collecting a subset of the added options
    mutable po::options_description* collectedOptions_;

    /// extracted command-line arguments are stored here
    po::variables_map variablesMap_;
};


// ----------------------------------------------------------------------------
// template definitions

template<class T>
void CommandLineParser::addOption(const std::string& name, OptionType type, const std::string& description) {
    if (type == RequiredOption)
        requiredOptions_.add_options()
#ifndef VRN_OLD_BOOST
            (name.c_str(), po::value<T>()->required(), description.c_str());
#else
            (name.c_str(), po::value<T>(), description.c_str());
#endif
    else if (type == MainOption)
        mainOptions_.add_options()
            (name.c_str(), po::value<T>(), description.c_str());
    else if (type == AdditionalOption)
        additionalOptions_.add_options()
            (name.c_str(), po::value<T>(), description.c_str());
    else
        tgtAssert(false, "unknown OptionType");
}

template<class T>
void CommandLineParser::addOption(const std::string& name, OptionType type, const std::string& description,
        const T& defaultValue, const std::string& defaultValueText) {

    if (!defaultValueText.empty()) {
        if (type == RequiredOption)
#ifndef VRN_OLD_BOOST
            requiredOptions_.add_options()
                (name.c_str(), po::value<T>()->default_value(defaultValue, defaultValueText)->required(), description.c_str());
#else
            requiredOptions_.add_options()
                (name.c_str(), po::value<T>()->default_value(defaultValue, defaultValueText), description.c_str());
#endif
        else if (type == MainOption)
            mainOptions_.add_options()
                (name.c_str(), po::value<T>()->default_value(defaultValue, defaultValueText), description.c_str());
        else if (type == AdditionalOption)
            additionalOptions_.add_options()
                (name.c_str(), po::value<T>()->default_value(defaultValue, defaultValueText), description.c_str());
        else
            tgtAssert(false, "unknown OptionType");
    }
    else {
        if (type == RequiredOption)
#ifndef VRN_OLD_BOOST
            requiredOptions_.add_options()
                (name.c_str(), po::value<T>()->default_value(defaultValue)->required(), description.c_str());
#else
            requiredOptions_.add_options()
                (name.c_str(), po::value<T>()->default_value(defaultValue), description.c_str());
#endif
        else if (type == MainOption)
            mainOptions_.add_options()
                (name.c_str(), po::value<T>()->default_value(defaultValue), description.c_str());
        else if (type == AdditionalOption)
            additionalOptions_.add_options()
                (name.c_str(), po::value<T>()->default_value(defaultValue), description.c_str());
        else
            tgtAssert(false, "unknown OptionType");
    }
}

template<class T>
void voreen::CommandLineParser::addOption(const std::string& name, T& out, OptionType type, const std::string& description) {
    if (type == RequiredOption)
#ifndef VRN_OLD_BOOST
        requiredOptions_.add_options()
            (name.c_str(), po::value<T>(&out)->required(), description.c_str());
#else
        requiredOptions_.add_options()
            (name.c_str(), po::value<T>(&out), description.c_str());
#endif
    else if (type == MainOption)
        mainOptions_.add_options()
            (name.c_str(), po::value<T>(&out), description.c_str());
    else if (type == AdditionalOption)
        additionalOptions_.add_options()
            (name.c_str(), po::value<T>(&out), description.c_str());
    else
        tgtAssert(false, "unknown OptionType");
}

template<class T>
void voreen::CommandLineParser::addOption(const std::string& name, T& out, OptionType type, const std::string& description,
        const T& defaultValue, const std::string& defaultValueText) {

    if (!defaultValueText.empty()) {
        if (type == RequiredOption)
#ifndef VRN_OLD_BOOST
            requiredOptions_.add_options()
                (name.c_str(), po::value<T>(&out)->default_value(defaultValue, defaultValueText)->required(), description.c_str());
#else
            requiredOptions_.add_options()
                (name.c_str(), po::value<T>(&out)->default_value(defaultValue, defaultValueText), description.c_str());
#endif
        else if (type == MainOption)
            mainOptions_.add_options()
                (name.c_str(), po::value<T>(&out)->default_value(defaultValue, defaultValueText), description.c_str());
        else if (type == AdditionalOption)
            additionalOptions_.add_options()
                (name.c_str(), po::value<T>(&out)->default_value(defaultValue, defaultValueText), description.c_str());
        else
            tgtAssert(false, "unknown OptionType");
    }
    else {
        if (type == RequiredOption)
#ifndef VRN_OLD_BOOST
            requiredOptions_.add_options()
                (name.c_str(), po::value<T>(&out)->default_value(defaultValue)->required(), description.c_str());
#else
            requiredOptions_.add_options()
                (name.c_str(), po::value<T>(&out)->default_value(defaultValue), description.c_str());
#endif
        else if (type == MainOption)
            mainOptions_.add_options()
                (name.c_str(), po::value<T>(&out)->default_value(defaultValue), description.c_str());
        else if (type == AdditionalOption)
            additionalOptions_.add_options()
                (name.c_str(), po::value<T>(&out)->default_value(defaultValue), description.c_str());
        else
            tgtAssert(false, "unknown OptionType");
    }
}

template<class T>
void CommandLineParser::addMultiOption(const std::string& name, OptionType type, const std::string& description) {
    if (type == RequiredOption)
#ifndef VRN_OLD_BOOST
        requiredOptions_.add_options()
            (name.c_str(), po::value<std::vector<T> >()->multitoken()->required(), description.c_str());
#else
        requiredOptions_.add_options()
            (name.c_str(), po::value<std::vector<T> >()->multitoken(), description.c_str());
#endif
    else if (type == MainOption)
        mainOptions_.add_options()
            (name.c_str(), po::value<std::vector<T> >()->multitoken(), description.c_str());
    else if (type == AdditionalOption)
        additionalOptions_.add_options()
            (name.c_str(), po::value<std::vector<T> >()->multitoken(), description.c_str());
    else
        tgtAssert(false, "unknown OptionType");
}

template<class T>
void CommandLineParser::addMultiOption(const std::string& name, std::vector<T>& out, OptionType type, const std::string& description) {
    if (type == RequiredOption)
#ifndef VRN_OLD_BOOST
        requiredOptions_.add_options()
            (name.c_str(), po::value<std::vector<T> >(&out)->multitoken()->required(), description.c_str());
#else
        requiredOptions_.add_options()
            (name.c_str(), po::value<std::vector<T> >(&out)->multitoken(), description.c_str());
#endif
    else if (type == MainOption)
        mainOptions_.add_options()
            (name.c_str(), po::value<std::vector<T> >(&out)->multitoken(), description.c_str());
    else if (type == AdditionalOption)
        additionalOptions_.add_options()
            (name.c_str(), po::value<std::vector<T> >(&out)->multitoken(), description.c_str());
    else
        tgtAssert(false, "unknown OptionType");
}

template<class T>
void CommandLineParser::getOptionValue(const std::string& name, T& out) const {
    if (isOptionSet(name))
        out = variablesMap_[name].as<T>();
}

template<class T>
void CommandLineParser::getMultiOptionValues(const std::string& name, std::vector<T>& out) const {
    if (isOptionSet(name))
        out = variablesMap_[name].as<std::vector<T> >();
}

} // namespace

#endif // VRN_COMMANDLINEPARSER_H
