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

#include "voreen/qt/aboutboxbase.h"

#include "voreen/qt/mainwindow/voreenqtmainwindow.h"

#include "voreen/core/version.h"

#include <QSizePolicy>
#include <QLabel>
#include <QVBoxLayout>
#include <QPushButton>

namespace voreen {

void AboutBoxBase::initSoftwareDescription() {
    softwareDescription_ = "Volume Rendering Engine " + VoreenVersion::getVersion() + "\n"
                           + VoreenVersion::getCopyright();
}

void AboutBoxBase::initHomepageRevisionText() {
    homepageRevisionText_ = "Revision: " + VoreenVersion::getRevision().substr(0,8) + "<br>"
                            + "<br>"
                            +"<a href=\"http://voreen.uni-muenster.de\">http://voreen.uni-muenster.de</a>";
}

void AboutBoxBase::initDevelopersList() {
    if(developers_.empty()) {
        developers_ << QString::fromLatin1("Emad Altamimi");
        developers_ << QString::fromLatin1("Alexander Bock");
        developers_ << QString::fromLatin1("Benjamin Bolte");
        developers_ << QString::fromLatin1("Helge Böschen");
        developers_ << QString::fromLatin1("Stephan Brandt");
        developers_ << QString::fromLatin1("Annika Bürger");
        developers_ << QString::fromLatin1("Raphael Bruns");
        developers_ << QString::fromLatin1("Mathias Dehne");
        developers_ << QString::fromLatin1("Christian Döring");
        developers_ << QString::fromLatin1("Matthias Droste");
        developers_ << QString::fromLatin1("Maike Dudek");
        developers_ << QString::fromLatin1("Maik Dworczynski");
        developers_ << QString::fromLatin1("Jan Esser");
        developers_ << QString::fromLatin1("André Exeler");
        developers_ << QString::fromLatin1("Björn Feischen");
        developers_ << QString::fromLatin1("Dirk Feldmann");
        developers_ << QString::fromLatin1("Alejandro Figueroa Meana");
        developers_ << QString::fromLatin1("Timo Griese");
        developers_ << QString::fromLatin1("Jeffrey Hall");
        developers_ << QString::fromLatin1("Philipp Hanraths");
        developers_ << QString::fromLatin1("Bernd Hemmer");
        developers_ << QString::fromLatin1("Dieter Janzen");
        developers_ << QString::fromLatin1("Jens Kasten");
        developers_ << QString::fromLatin1("Daniel Kirsch");
        developers_ << QString::fromLatin1("Florian Kleene");
        developers_ << QString::fromLatin1("Benjamin König");
        developers_ << QString::fromLatin1("Rico Lehmann");
        developers_ << QString::fromLatin1("Roland Leißa");
        developers_ << QString::fromLatin1("Simon Leistikow");
        developers_ << QString::fromLatin1("Sören Linnemann");
        developers_ << QString::fromLatin1("Markus Madeja");
        developers_ << QString::fromLatin1("Zoha Moztarzadeh");
        developers_ << QString::fromLatin1("Reza Nawrozi");
        developers_ << QString::fromLatin1("Borislav Petkov");
        developers_ << QString::fromLatin1("Carsten Praßni");
        developers_ << QString::fromLatin1("Stephan Rademacher");
        developers_ << QString::fromLatin1("Eelamayooran Raveendran");
        developers_ << QString::fromLatin1("Rainer Reich");
        developers_ << QString::fromLatin1("Mona Riemenschneider");
        developers_ << QString::fromLatin1("Christoph Rosemann");
        developers_ << QString::fromLatin1("Jan Roters");
        developers_ << QString::fromLatin1("Sönke Schmid");
        developers_ << QString::fromLatin1("Christian Schulte zu Berge");
        developers_ << QString::fromLatin1("Yannik Siegert");
        developers_ << QString::fromLatin1("Michael Specht");
        developers_ << QString::fromLatin1("Fabian Spiegel");
        developers_ << QString::fromLatin1("Sven Strothoff");
        developers_ << QString::fromLatin1("Tahar Talebi");
        developers_ << QString::fromLatin1("Sebastian Terhorst");
        developers_ << QString::fromLatin1("David Terbeek");
        developers_ << QString::fromLatin1("Alexander Theißen");
        developers_ << QString::fromLatin1("Nils Vensler");
        developers_ << QString::fromLatin1("Andreas Völker");
        developers_ << QString::fromLatin1("Christian Vorholt");
        developers_ << QString::fromLatin1("Carolin Walter");
        developers_ << QString::fromLatin1("Paul Weingardt");
        developers_ << QString::fromLatin1("Michael Weinkath");
        developers_ << QString::fromLatin1("Sascha Wendt");
        developers_ << QString::fromLatin1("Malte Wildt");
        developers_ << QString::fromLatin1("Frank Wisniewski");
        developers_ << QString::fromLatin1("Marco Ziolkowski");
        developers_ << QString::fromLatin1("Johannes Zurhorst");
    }
}

void AboutBoxBase::initMainDevelopersList() {
    if(mainDevelopers_.empty()) {
        mainDevelopers_ << QString::fromLatin1("Tobias Brix");
        mainDevelopers_ << QString::fromLatin1("Stefan Diepenbrock");
        mainDevelopers_ << QString::fromLatin1("Dominik Drees");
        mainDevelopers_ << QString::fromLatin1("Florian Lindemann");
        mainDevelopers_ << QString::fromLatin1("Jörg Mensmann");
        mainDevelopers_ << QString::fromLatin1("Jennis Meyer-Spradow");
        mainDevelopers_ << QString::fromLatin1("Jörg-Stefan Praßni");
        mainDevelopers_ << QString::fromLatin1("Timo Ropinski");
        mainDevelopers_ << QString::fromLatin1("Aaron Scherzinger");
    }
}

void AboutBoxBase::initLicenseText() {
    licenseString_ = "You may use, distribute and copy the Voreen software package "\
                     "under the terms of the GNU General Public License version 2, "\
                     "see the files LICENSE.txt and LICENSE-academic.txt for details.";
}


AboutBoxBase::AboutBoxBase(VoreenQtMainWindow* mainWindow, QString imagePath)
    : QDialog(mainWindow)
    , initialized_(false)
    , mainWindow_(mainWindow)
    , imagePath_(imagePath)
{
    setWindowFlags(windowFlags() | Qt::MSWindowsFixedSizeDialogHint); //removes resize mouse event on hover
    //set color background
#if (QT_VERSION >= 0x040400) && !defined(__APPLE__) && !defined(VRN_NO_STYLESHEET)
    setStyleSheet("QDialog { background-color: qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 #444444, stop:1 #aaaaaa) }");
#endif
}

void AboutBoxBase::initialize() {
    if(initialized_) return;

    initSoftwareDescription();
    initHomepageRevisionText();
    initDevelopersList();
    initMainDevelopersList();
    initLicenseText();

    initAndLayoutItems();

    initialized_ = true;
}

int AboutBoxBase::exec() {
    if(!initialized_)
        initialize();
    return QDialog::exec();
}

void AboutBoxBase::initAndLayoutItems() {
    //used font
    QFont baseFont(QString("Sans Serif"),9,400);

    //set icon and title
    setWindowIcon(mainWindow_->windowIcon());
    setWindowTitle(QString(mainWindow_->getApplicationTitle().c_str()));

    QVBoxLayout* mainLayout = new QVBoxLayout(this);

    //set image
    QPixmap image =  QPixmap(imagePath_);
    QLabel* imageLabel = new QLabel(this);
    imageLabel->setMinimumSize(QSize(510, 86));
    imageLabel->setLineWidth(0);
    imageLabel->setPixmap(QPixmap(imagePath_));
    imageLabel->setAlignment(Qt::AlignCenter);
    mainLayout->addWidget(imageLabel);

    //set top information
    QHBoxLayout* topLayout = new QHBoxLayout();
    QLabel* softwareDescriptionLabel = new QLabel(this);
    softwareDescriptionLabel->setAlignment(static_cast<Qt::Alignment>(Qt::AlignLeft|Qt::AlignTop));
    softwareDescriptionLabel->setText(QString::fromLatin1(softwareDescription_.c_str()));
    QLabel* homepageRevisionLabel = new QLabel(this);
    homepageRevisionLabel->setFont(baseFont);
    homepageRevisionLabel->setAlignment(static_cast<Qt::Alignment>(Qt::AlignRight|Qt::AlignTop));
    homepageRevisionLabel->setOpenExternalLinks(true);
    homepageRevisionLabel->setText(tr("<font color=\"black\">") + QString(homepageRevisionText_.c_str()) + tr("<\\font>"));

    topLayout->addWidget(softwareDescriptionLabel);
    topLayout->addStretch();
    topLayout->addWidget(homepageRevisionLabel);

    mainLayout->addLayout(topLayout);

    //add developers
    QLabel* mainDevelopersLabel = new QLabel(QString(std::string("<span style=\" font-weight:600;\">Concept & Design: </span>"
                                                             + std::string("<font color=\"black\">") + convertStringListToString(mainDevelopers_) + std::string("<\\font>")).c_str()),this);
    mainDevelopersLabel->setWordWrap(true);
    mainDevelopersLabel->setFont(baseFont);
    mainLayout->addWidget(mainDevelopersLabel);

    //add developers
    QLabel* developersLabel = new QLabel(QString(std::string("<span style=\" font-weight:600;\">Developers: </span>"
                                                             + std::string("<font color=\"black\">") + convertStringListToString(developers_) + std::string("<\\font>")).c_str()),this);
    developersLabel->setWordWrap(true);
    developersLabel->setFont(baseFont);
    mainLayout->addWidget(developersLabel);

    //add licence
    QLabel* licenseLabel = new QLabel(QString(std::string("<span style=\" font-weight:600;\">License: </span>" + std::string("<font color=\"black\">") + licenseString_ + std::string("<\\font>")).c_str()),this);
    licenseLabel->setWordWrap(true);
    licenseLabel->setFont(baseFont);
    mainLayout->addWidget(licenseLabel);

    // add button
    QHBoxLayout* buttonLayout = new QHBoxLayout();
    buttonLayout->addStretch();
    QPushButton* okButton = new QPushButton("OK",this);
    okButton->setDefault(true);
    buttonLayout->addWidget(okButton);
    buttonLayout->addStretch();
    mainLayout->addLayout(buttonLayout);

    // size policy
    setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Fixed);

    layout()->setSizeConstraint( QLayout::SetFixedSize );
    setSizeGripEnabled(false);

    adjustSize();

    QObject::connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
}

std::string AboutBoxBase::convertStringListToString(const QStringList& list) {
    std::string result("");
    for (int i=0; i < list.size(); i++) {
        if (i > 0)
            result += ", ";
        result += list[i].toStdString();
    }
    return result;
}

} // namespace
