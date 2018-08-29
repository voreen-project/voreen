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

#ifndef VRN_TRANSFUNCPLUGINBASE_H
#define VRN_TRANSFUNCPLUGINBASE_H

#include <QWidget>

namespace voreen {

    class TransFuncProperty;
    class TransFuncEditorBase;
    class TransFuncWidgetBase;

/**
 * Container class for transfer function editors and widgets. Each transfer function has to
 * implement this plugin class. The transfer function property widget only used this interface to communicate.
 */
class TransFuncPluginBase : public QWidget {
    Q_OBJECT

        friend class TransFuncPropertyWidget;
public:

    /** Constructor */
    TransFuncPluginBase(TransFuncProperty* prop, QWidget* parent = 0);

    /** Destructor */
    ~TransFuncPluginBase();

    /** has to be implemented in the sub-class. */
    virtual TransFuncEditorBase* createEditor() = 0;
    /** has to be implemented in the sub-class. */
    virtual TransFuncWidgetBase* createWidget() = 0;

protected:


private:
    TransFuncProperty* property_;    ///< transfer function property that belongs to this plugin

    TransFuncEditorBase* editor_;    ///< editor used in this plugin
    TransFuncWidgetBase* widget_;    ///< widget used in this plugin
};

} // namespace voreen

#endif // VRN_TRANSFUNCPLUGINBASE_H
