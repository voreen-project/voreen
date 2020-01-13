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

#include "voreen/qt/networkeditor/meta/textboxmetadata.h"

#include "voreen/qt/networkeditor/graphicitems/textboxes/textboxbasegraphicsitem.h"

#include "voreen/core/utils/stringutils.h"

namespace voreen {

TextBoxMetaData::TextBoxMetaData(const NWEBaseGraphicsItemUserTypes userType, const std::string caption, const int captionFontSize,
                                 const std::string content, const int contentFontSize, const int x, const int y, const int width, const int height,
                                 const QColor fontColor, const QColor baseColor, const bool showCaption)
    : userType_(userType), caption_(caption), captionFontSize_(captionFontSize), content_(content), contentFontSize_(contentFontSize)
    , x_(x), y_(y), width_(width), height_(height), fontColor_(fontColor), baseColor_(baseColor), showCaption_(showCaption)
{
}

TextBoxMetaData::TextBoxMetaData(const TextBoxBaseGraphicsItem* item)
{
    tgtAssert(item, "no item passed");
    userType_ = static_cast<NWEBaseGraphicsItemUserTypes>(item->type());
    caption_ = "";
    captionFontSize_ = -1;
    if(const QGraphicsTextItem* ti = item->getCaptionItem()) {
        caption_ = ti->toPlainText().toStdString();
        captionFontSize_ = ti->font().pointSize();
    }
    content_ = "";
    contentFontSize_ = -1;
    if(const QTextEdit* te = item->getContentEditor()) {
        content_ = te->toPlainText().toStdString();
        contentFontSize_ = te->font().pointSize();
    }
    x_ = item->pos().x();
    y_ = item->pos().y();
    width_ = item->getCurrentSize().width();
    height_ = item->getCurrentSize().height();
    fontColor_ = item->getFontColor();
    baseColor_ = item->getBaseColor();
    showCaption_ = item->isCaptionEnabled();
}

std::string TextBoxMetaData::toString() const {
    std::stringstream s;
    s << "type = " << itos(userType_) << "; ";

    s << "caption = " << caption_ << ", ";
    s << "captionFontSize = " << captionFontSize_ << ", ";

    s << "content = " << content_ << "; ";
    s << "contentFontSize = " << contentFontSize_ << ", ";

    s << "x = " << itos(x_) << ", ";
    s << "y = " << itos(y_) << "; ";

    s << "width = " << itos(width_) << ", ";
    s << "height = " << itos(height_) << "; ";

    s << "fontColor = " << fontColor_.name().toStdString() << ", ";
    s << "baseColor = " << baseColor_.name().toStdString() << "; ";

    s << "showCaption = " << (showCaption_ ? "true" : "false") << ";";

    return s.str();
}

std::string TextBoxMetaData::toString(const std::string& component) const {
    if (component == "type")
        return itos(userType_);
    else if (component == "caption")
        return caption_;
    else if (component == "captionFontSize")
        return itos(captionFontSize_);
    else if (component == "content")
        return content_;
    else if (component == "contentFontSize")
        return itos(contentFontSize_);
    else if (component == "x")
        return itos(x_);
    else if (component == "y")
        return itos(y_);
    else if (component == "width")
        return itos(width_);
    else if (component == "height")
        return itos(height_);
    else if (component == "fontColor")
        return fontColor_.name().toStdString();
    else if (component == "baseColor")
        return baseColor_.name().toStdString();
    else if (component == "showCaption")
        return (showCaption_ ? "true" : "false");
    else return toString();
}

void TextBoxMetaData::serialize(Serializer& s) const {
    s.serialize("userType", static_cast<int>(userType_));
    if(!caption_.empty()) s.serialize("caption",caption_);
    if(captionFontSize_ != -1) s.serialize("captionFontSize",captionFontSize_);
    if(!content_.empty()) s.serialize("content",content_);
    if(contentFontSize_ != -1) s.serialize("contentFontSize",contentFontSize_);
    if(x_ != -1) s.serialize("x", x_);
    if(y_ != -1) s.serialize("y", y_);
    if(width_ != -1) s.serialize("width", width_);
    if(height_ != -1) s.serialize("height", height_);
    if(fontColor_ != Qt::black) s.serialize("fontColor", fontColor_.name().toStdString());
    if(baseColor_ != Qt::gray) s.serialize("baseColor", baseColor_.name().toStdString());
    if(!showCaption_) s.serialize("showCaption",showCaption_);
}

void TextBoxMetaData::deserialize(Deserializer& s) {
    int tmp;
    s.deserialize("userType", tmp);
    userType_ = static_cast<NWEBaseGraphicsItemUserTypes>(tmp);
    try {
        s.deserialize("caption", caption_);
    } catch (SerializationNoSuchDataException&) {
        // caption was not serialized so just ignore...
        s.removeLastError();
    }
    try {
        s.deserialize("captionFontSize", captionFontSize_);
    } catch (SerializationNoSuchDataException&) {
        // captionFontSize was not serialized so just ignore...
        s.removeLastError();
        captionFontSize_ = 9;
    }
    try {
        s.deserialize("content", content_);
    } catch (SerializationNoSuchDataException&) {
        // content was not serialized so just ignore...
        s.removeLastError();
    }
    try {
        s.deserialize("contentFontSize", contentFontSize_);
    } catch (SerializationNoSuchDataException&) {
        // contentFontSize was not serialized so just ignore...
        s.removeLastError();
        contentFontSize_ = 6;
    }
    try {
        s.deserialize("x", x_);
    } catch (SerializationNoSuchDataException&) {
        // x position was not serialized so just ignore...
        s.removeLastError();
    }
    try {
        s.deserialize("y", y_);
    } catch (SerializationNoSuchDataException&) {
        // y position was not serialized so just ignore...
        s.removeLastError();
    }
    try {
        s.deserialize("width", width_);
    } catch (SerializationNoSuchDataException&) {
        // width was not serialized so just ignore...
        s.removeLastError();
    }
    try {
        s.deserialize("height", height_);
    } catch (SerializationNoSuchDataException&) {
        // height was not serialized so just ignore...
        s.removeLastError();
    }
    try {
        std::string tmp;
        s.deserialize("fontColor", tmp);
        fontColor_ = QColor(tmp.c_str());
    } catch (SerializationNoSuchDataException&) {
        // font color was not serialized so just ignore...
        s.removeLastError();
    }
    try {
        std::string tmp;
        s.deserialize("baseColor", tmp);
        baseColor_ = QColor(tmp.c_str());
    } catch (SerializationNoSuchDataException&) {
        // base color was not serialized so just ignore...
        s.removeLastError();
    }
     try {
        s.deserialize("showCaption", showCaption_);
    } catch (SerializationNoSuchDataException&) {
        // show caption was not serialized so just ignore...
        s.removeLastError();
    }
}

void TextBoxMetaData::setUserType(const NWEBaseGraphicsItemUserTypes userType) {
    userType_ = userType;
}
NWEBaseGraphicsItemUserTypes TextBoxMetaData::getUserType() const {
    return userType_;
}


void TextBoxMetaData::setCaption(const std::string caption) {
    caption_ = caption;
}
std::string TextBoxMetaData::getCaption() const {
    return caption_;
}


void TextBoxMetaData::setCaptionFontSize(const int fontSize) {
    captionFontSize_ = fontSize;
}
int TextBoxMetaData::getCaptionFontSize() const {
    return captionFontSize_;
}

void TextBoxMetaData::setContent(const std::string content) {
    content_ = content;
}
std::string TextBoxMetaData::getContent() const {
    return content_;
}


void TextBoxMetaData::setContentFontSize(const int fontSize) {
    contentFontSize_ = fontSize;
}
int TextBoxMetaData::getContentFontSize() const {
    return contentFontSize_;
}

void TextBoxMetaData::setX(const int value) {
    x_ = value;
}
int TextBoxMetaData::getX() const {
    return x_;
}


void TextBoxMetaData::setY(const int value) {
    y_ = value;
}
int TextBoxMetaData::getY() const {
    return y_;
}


void TextBoxMetaData::setWidth(const int value) {
    width_ = value;
}
int TextBoxMetaData::getWidth() const {
    return width_;
}


void TextBoxMetaData::setHeight(const int value) {
    height_ = value;
}
int TextBoxMetaData::getHeight() const {
    return height_;
}

void TextBoxMetaData::setFontColor(const QColor col) {
    fontColor_ = col;
}
QColor TextBoxMetaData::getFontColor() const {
    return fontColor_;
}

void TextBoxMetaData::setBaseColor(const QColor col) {
    baseColor_ = col;
}
QColor TextBoxMetaData::getBaseColor() const {
    return baseColor_;
}

void TextBoxMetaData::setShowCaption(const bool b) {
    showCaption_ = b;
}
bool TextBoxMetaData::getShowCaption() const {
    return showCaption_;
}

MetaDataBase* TextBoxMetaData::clone() const {
    return new TextBoxMetaData(userType_, caption_, captionFontSize_, content_, contentFontSize_, x_, y_, width_, height_, fontColor_, baseColor_, showCaption_);
}

MetaDataBase* TextBoxMetaData::create() const {
    return new TextBoxMetaData();
}

} // namespace
