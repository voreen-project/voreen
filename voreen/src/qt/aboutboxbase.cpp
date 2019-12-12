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
        developers_ << "Emad Altamimi";
        developers_ << "Alexander Bock";
        developers_ << "Benjamin Bolte";
        developers_ << "Helge Böschen";
        developers_ << "Stephan Brandt";
        developers_ << "Annika Bürger";
        developers_ << "Raphael Bruns";
        developers_ << "Mathias Dehne";
        developers_ << "Christian Döring";
        developers_ << "Matthias Droste";
        developers_ << "Maike Dudek";
        developers_ << "Maik Dworczynski";
        developers_ << "Jan Esser";
        developers_ << "André Exeler";
        developers_ << "Björn Feischen";
        developers_ << "Dirk Feldmann";
        developers_ << "Alejandro Figueroa Meana";
        developers_ << "Timo Griese";
        developers_ << "Jeffrey Hall";
        developers_ << "Philipp Hanraths";
        developers_ << "Bernd Hemmer";
        developers_ << "Dieter Janzen";
        developers_ << "Jens Kasten";
        developers_ << "Daniel Kirsch";
        developers_ << "Florian Kleene";
        developers_ << "Benjamin König";
        developers_ << "Rico Lehmann";
        developers_ << "Roland Leißa";
        developers_ << "Sören Linnemann";
        developers_ << "Markus Madeja";
        developers_ << "Zoha Moztarzadeh";
        developers_ << "Reza Nawrozi";
        developers_ << "Borislav Petkov";
        developers_ << "Carsten Praßni";
        developers_ << "Stephan Rademacher";
        developers_ << "Eelamayooran Raveendran";
        developers_ << "Rainer Reich";
        developers_ << "Mona Riemenschneider";
        developers_ << "Christoph Rosemann";
        developers_ << "Jan Roters";
        developers_ << "Sönke Schmid";
        developers_ << "Christian Schulte zu Berge";
        developers_ << "Yannik Siegert";
        developers_ << "Michael Specht";
        developers_ << "Fabian Spiegel";
        developers_ << "Sven Strothoff";
        developers_ << "Tahar Talebi";
        developers_ << "Sebastian Terhorst";
        developers_ << "David Terbeek";
        developers_ << "Alexander Theißen";
        developers_ << "Nils Vensler";
        developers_ << "Andreas Völker";
        developers_ << "Christian Vorholt";
        developers_ << "Carolin Walter";
        developers_ << "Paul Weingardt";
        developers_ << "Michael Weinkath";
        developers_ << "Sascha Wendt";
        developers_ << "Malte Wildt";
        developers_ << "Frank Wisniewski";
        developers_ << "Marco Ziolkowski";
        developers_ << "Johannes Zurhorst";
    }
}

void AboutBoxBase::initMainDevelopersList() {
    if(mainDevelopers_.empty()) {
        mainDevelopers_ << "Tobias Brix";
        mainDevelopers_ << "Stefan Diepenbrock";
        mainDevelopers_ << "Dominik Drees";
        mainDevelopers_ << "Marina Evers";
        mainDevelopers_ << "Simon Leistikow";
        mainDevelopers_ << "Florian Lindemann";
        mainDevelopers_ << "Jörg Mensmann";
        mainDevelopers_ << "Jennis Meyer-Spradow";
        mainDevelopers_ << "Jörg-Stefan Praßni";
        mainDevelopers_ << "Timo Ropinski";
        mainDevelopers_ << "Aaron Scherzinger";
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
    softwareDescriptionLabel->setText(softwareDescription_.c_str());
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
