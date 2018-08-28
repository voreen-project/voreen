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

#ifndef VRN_TEXTBOXMETADATA_H
#define VRN_TEXTBOXMETADATA_H

#include "voreen/core/datastructures/meta/metadatabase.h"

#include "voreen/core/io/serialization/xmlserializer.h"
#include "voreen/core/io/serialization/xmldeserializer.h"

#include "voreen/qt/networkeditor/editor_settings.h"

#include "voreen/qt/voreenqtapi.h"

#include <QColor>

namespace voreen {

    class TextBoxBaseGraphicsItem;

/**
 * The @c TextBoxMetaData class stores all information of textboxgraphicsitems.
 *
 * @see MetaDataBase
 */
class VRN_QT_API TextBoxMetaData : public MetaDataBase {
public:
    /**
     * Creates a @c TextBoxMetaData object storing the given text graphics item.
     */
    TextBoxMetaData(const NWEBaseGraphicsItemUserTypes userType = UserTypesTextBoxBaseGraphicsItem, const std::string caption = "", const int captionFontSize = -1,
                    const std::string content = "", const int contentFontSize = -1, const int x = -1, const int y = -1, const int width = -1, const int height = -1,
                    const QColor fontColor = Qt::black, const QColor baseColor = Qt::gray, const bool showCaption = true);
    TextBoxMetaData(const TextBoxBaseGraphicsItem* item);
    virtual ~TextBoxMetaData() {}

    virtual std::string getClassName() const { return "TextBoxMetaData"; }
    virtual MetaDataBase* create() const;
    virtual MetaDataBase* clone() const;
    virtual std::string toString() const;
    virtual std::string toString(const std::string& /*component*/) const;

    /** @see Serializable::serialize */
    virtual void serialize(Serializer& s) const;
    /** @see Serializable::deserialize */
    virtual void deserialize(Deserializer& s);

    //----------------------------------------
    //  getter and setter
    //----------------------------------------

    void setUserType(const NWEBaseGraphicsItemUserTypes type);
    NWEBaseGraphicsItemUserTypes getUserType() const;

    void setCaption(const std::string caption);
    std::string getCaption() const;

    void setCaptionFontSize(const int fontSize);
    int getCaptionFontSize() const;

    void setContent(const std::string content);
    std::string getContent() const;

    void setContentFontSize(const int fontSize);
    int getContentFontSize() const;

    void setX(const int value);
    int getX() const;

    void setY(const int value);
    int getY() const;

    void setWidth(const int value);
    int getWidth() const;

    void setHeight(const int value);
    int getHeight() const;

    void setFontColor(const QColor col);
    QColor getFontColor() const;

    void setBaseColor(const QColor col);
    QColor getBaseColor() const;

    void setShowCaption(const bool b);
    bool getShowCaption() const;

private:
    NWEBaseGraphicsItemUserTypes userType_;
    std::string caption_;
    int captionFontSize_;
    std::string content_;
    int contentFontSize_;
    int x_;
    int y_;
    int width_;
    int height_;
    QColor fontColor_;
    QColor baseColor_;
    bool showCaption_;
};

} // namespace

#endif // VRN_TEXTBOXMETADATA_H
