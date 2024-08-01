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

#include "voreen/core/utils/commandlineparser.h"

#include "tgt/filesystem.h"
#include "voreen/core/utils/stringutils.h"

#include <vector>
#include <map>
#include <stdlib.h>

using std::vector;
using std::pair;
using std::string;

namespace voreen {

CommandLineParser::CommandLineParser(const std::string& programName, const std::string& programDescription,
        po::command_line_style::style_t cmdStyle, const std::string& configFile)
  : executableName_(programName)
  , programDescription_(programDescription)
  , cmdStyle_(cmdStyle)
  , configFile_(configFile)
  , requiredOptions_("Required Options")
  , mainOptions_("Main Options")
  , additionalOptions_("Additional Options")
  , collectedOptions_(0)
  , variablesMap_()
{}

CommandLineParser::~CommandLineParser() {
    delete collectedOptions_;
    collectedOptions_ = 0;
}

void CommandLineParser::setCommandLine(int argc, char** argv) {
    // argv[0] = program name
    // argv[i] = i-th argument
    if (argc > 0 && argv && argv[0])
        programPath_ = argv[0];
    else
        programPath_ = "";

    argc_ = argc;
    argv_ = argv;
}

void CommandLineParser::setCommandLineStyle(po::command_line_style::style_t cmdStyle) {
    cmdStyle_ = cmdStyle;
}

void CommandLineParser::setConfigFile(const std::string& configFile) {
    configFile_ = configFile;
}

void CommandLineParser::execute() {

    const po::options_description* optionCollection = getOptions(AllTypes);
    try {
        // parse command line
        po::store(po::parse_command_line(argc_, argv_, *optionCollection, cmdStyle_), variablesMap_);

        // parse config file (options extracted from command line are not overridden)
        if (!configFile_.empty() && tgt::FileSystem::fileExists(configFile_)) {
#ifndef VRN_OLD_BOOST
            std::cout << "Reading config file: " << tgt::FileSystem::cleanupPath(configFile_) << std::endl;
            po::store(po::parse_config_file<char>(configFile_.c_str(), *optionCollection, true), variablesMap_);
#endif
        }

        po::notify(variablesMap_);
    }
    catch (std::exception& e) {
        throw VoreenException(e.what());
    }
}


void CommandLineParser::addFlagOption(const std::string& name, OptionType type, const std::string& description) {
    if (type == RequiredOption)
#ifndef VRN_OLD_BOOST
        requiredOptions_.add_options()
            (name.c_str(), po::value<bool>()->zero_tokens()->required(), description.c_str());
#else
        requiredOptions_.add_options()
            (name.c_str(), po::value<bool>()->zero_tokens(), description.c_str());
#endif
    else if (type == MainOption)
        mainOptions_.add_options()
            (name.c_str(), po::value<bool>()->zero_tokens(), description.c_str());
    else if (type == AdditionalOption)
        additionalOptions_.add_options()
            (name.c_str(), po::value<bool>()->zero_tokens(), description.c_str());
}

void CommandLineParser::addFlagOption(const std::string& name, bool& out, OptionType type, const std::string& description) {
    if (type == RequiredOption)
#ifndef VRN_OLD_BOOST
        requiredOptions_.add_options()
            (name.c_str(), po::value<bool>(&out)->zero_tokens()->required(), description.c_str());
#else
        requiredOptions_.add_options()
            (name.c_str(), po::value<bool>(&out)->zero_tokens(), description.c_str());
#endif
    else if (type == MainOption)
        mainOptions_.add_options()
            (name.c_str(), po::value<bool>(&out)->zero_tokens(), description.c_str());
    else if (type == AdditionalOption)
        additionalOptions_.add_options()
            (name.c_str(), po::value<bool>(&out)->zero_tokens(), description.c_str());
}

std::string CommandLineParser::getUsageString(OptionType optionTypes, bool programDescription) const {
    std::ostringstream stream;
    if (programDescription && !programDescription_.empty())
        stream << programDescription_ << "\n\n";
    stream << "Usage: " << executableName_ << " [options]" << "\n";
    stream << *getOptions(optionTypes);
    return stream.str();
}

const boost::program_options::options_description* CommandLineParser::getOptions(OptionType types) const {
    delete collectedOptions_;
    collectedOptions_ = new po::options_description();

    if ((types & RequiredOption) && !requiredOptions_.options().empty())
        collectedOptions_->add(requiredOptions_);
    if ((types & MainOption) && !mainOptions_.options().empty())
        collectedOptions_->add(mainOptions_);
    if ((types & AdditionalOption) && !additionalOptions_.options().empty())
        collectedOptions_->add(additionalOptions_);

    return collectedOptions_;
}

bool CommandLineParser::isOptionSet(const std::string& name) const {
    return variablesMap_.count(name);
}

std::string CommandLineParser::getProgramPath() const {
    return programPath_;
}

} // namespace
