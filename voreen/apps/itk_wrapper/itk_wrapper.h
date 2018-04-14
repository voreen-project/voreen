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

#ifndef ITK_WRAPPER_H
#define ITK_WRAPPER_H

#include "baseclasses.h"

namespace voreen {

/**
 * An ITK_Wrapper provides methods to generate Voreen-processors representing ITK-filters.
*/
class ITK_Wrapper{
public:
    ITK_Wrapper();
    ITK_Wrapper(std::vector<std::string> modules);

    /**
     * Initializes some string-vectors (e.g. codeStates_) with possible values
     * for filter variables (e.g. codeState).
     * (These are used to except errors in the xml-files of the filtermodules.)
     */
    void initializeConstants();

    /**
     * Deserializes the xml-files of the filtermodules and checks each filter for correct
     * variable-setting. If no error occurred, itk_generatedmodule.cpp,
     * itk_generated_module.xml and all filter-.cpp and -.h files will be generated.
     */
    void wrapITKtoVoreen();

    static const std::string loggerCat_;


private:
    /**
     * Deserializes the xml-files of the filtermodules and checks each module/filter for correct
     * variable-setting by calling the checkFilters-method.
     */
    void deserializeModuleXMLs();

    /**
     * Checks each filter in module itkModule for correct variable-setting and produces
     * an error-massage if a variable has a value which is not supported or makes no sense.
     * (In this way errors in the xml-files of the filtermodules can be found
     * and eliminated.)
     *
     * @param itkModule, module to be checked
     */
    void checkFilters(ITK_Module itkModule);

    /**
     * Creates itk_generatedmodule.cpp
     */
    void createITKGeneratedModule_CPP();

    /**
     * Creates itk_generated.cmake
     */
    void createITKGenerated_CMake();

    /**
     * Creates all filter .cpp- and .h-files
     */
    void createFilters();

    /**
     * Returns the path of the directory containing the ITK-xml-files
     *
     * @return path, path of the directory containing the ITK-xml-files
     */
    std::string getXMLPath();

    /**
     * Returns the path of a special ITK-xml-file
     *
     * @param xmlFilename
     * @return path, path of the ITK-xml-file with the filname xmlFilename
     */
    std::string getXMLPath(std::string xmlFilename);

    /**
     * Returns the path of the ITK-module-directories
     *
     * @return path, path of the ITK-module-directories
     */
    std::string getProcessorOutputPath();

    /**
     * Returns the path of a special ITK-module-directory where the
     * processors of the module should be stored.
     *
     * @param moduleName
     * @return path, path of the ITK-module-directory with the name moduleName
     */
    std::string getProcessorOutputPath(std::string moduleName);

    /**
     * Returns the path of a special template-file.
     *
     * @param templateFilename
     * @return path, path of the template with the filename templateFilename
     */
    std::string getTemplatePath(std::string templateFilename);

    /**
     * Returns the path of a special filter in the specialFilters directory.
     *
     * @param specialFilterFilename
     * @return path, path of the special filter with the filename specialFilterFilename
     */
    std::string getSpecialFilterPath(std::string specialFilterFilename);

    std::vector<std::string> moduleXMLs_;           ///< Contains xml-filenames of modules to be wrapped
    std::vector<std::string> codeStates_;           ///< Contains possible "values" for codeState_ of a filter
    std::vector<std::string> integerScalarTypes_;   ///< Contains possible scalar-integer volumtypes for possibleTypes_ of a port
    std::vector<std::string> realScalarTypes_;      ///< Contains possible scalar-real volumtypes for possibleTypes_ of a port
    std::vector<std::string> scalarTypes_;          ///< Contains possible scalar volumtypes for possibleTypes_ of a port
    std::vector<std::string> integerVectorTypes_;   ///< Contains possible vector-integer volumtypes for possibleTypes_ of a port
    std::vector<std::string> realVectorTypes_;      ///< Contains possible vector-real volumtypes for possibleTypes_ of a port
    std::vector<std::string> vectorTypes_;          ///< Contains possible vector volumtypes for possibleTypes_ of a port
    std::vector<std::string> allTypes_;             ///< Contains all possible volumtypes for possibleTypes_ of a port
    std::vector<std::string> metaTypes_;            ///< Contains possible meta-types for possibleTypes_ of a port
    std::vector<std::string> argumentTypes_;        ///< Contains possible types for type_ of an argument
    std::vector<ITK_Module> filtermodules_;         ///< Contains all filtermodules with all filters to be generated
};

} //namespace voreen

#endif // ITK_WRAPPER_H
