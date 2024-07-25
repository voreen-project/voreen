/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Simon Großmann
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef CSV_WRITER_HH
#define CSV_WRITER_HH

#include <fstream>
#include <iostream>
#include <unistd.h>
#include "core/singleton.h"
#include "io/fileName.h"
#include "CSVWriter.h"
#include "utilities/vectorHelpers.h"

namespace olb {
/// Constructor with name of outputFiles
/// boolean true for real-time plotting //WARNING: experimental!
template< typename T >
CSV<T>::CSV(std::string name, char separator, std::vector<std::string> columnTags)
  : _name(name),
    _dataFile(singleton::directories().getGnuplotOutDir()+"data/"+_name+".dat"),
    _dir(singleton::directories().getGnuplotOutDir()),
    _separator(separator),
    _columnTags(columnTags)
{
  if (singleton::mpi().getRank() == _rank) {
    std::ofstream fout;

    ///add (new) data file
    fout.open(_dataFile.c_str(), std::ios::trunc);
    if(columnTags.size()>0)
    {
      setColumnTags(columnTags);
    }
    fout.close();
  }
}

template< typename T >
CSV<T>::CSV(std::string name, std::vector<std::string> columnTags) : CSV(name, ' ', columnTags){}

template< typename T >
CSV<T>::CSV(std::string name) : CSV(name, ' ', std::vector<std::string> {}){}

template< typename T >
CSV<T>::CSV(std::string name, char separator) : CSV(name, separator, std::vector<std::string> {}){}

/// writes the data file for one double and a vector of doubles with a specific filename
template< typename T >
void CSV<T>::writeDataFile(T xValue, const std::vector<T>& yValues,const std::string& plotNameFile, int precision)
{
  if (singleton::mpi().getRank() == _rank) {
    std::ofstream fout;
    std::string DATAF;
    DATAF = singleton::directories().getGnuplotOutDir()+"data/"+plotNameFile+".dat";
    fout.precision(precision);
    fout.open(DATAF.c_str(), std::ios::out | std::ios::app);
    fout << xValue;
    for (unsigned int i = 0; i < yValues.size(); i++) {
      fout << _separator << yValues[i];
    }
    fout << std::endl;
    fout.close();
  }
  return;
}


/// writes the data file for two doubles with a specific filename
template< typename T >
void CSV<T>::writeDataFile(T xValue, T yValue, const std::string& plotNameFile, int precision)
{
  std::vector<T> yValues{yValue};
  writeDataFile(xValue, yValues, plotNameFile, precision);
  return;
}


/// writes the data file for one double and a vector of doubles (x and y1,y2,...)
template< typename T >
void CSV<T>::writeDataFile(T xValue, const std::vector<T>& yValues, int precision)
{
  writeDataFile(xValue, yValues, _name, precision);
  return;
}

/// writes the data file for two doubles (x and y)
template< typename T >
void CSV<T>::writeDataFile(T xValue, T yValue, int precision)
{
  std::vector<T> yValues{yValue};
  writeDataFile(xValue, yValues, _name,precision);
  return;
}

/// adds column tags at the beginning of the csv data file.
/// saves the content of the file and rewrites it with the tags at the top and then inserting the content again
template< typename T >
void CSV<T>::setColumnTags(const std::vector<std::string> columnTags, std::string& plotFileName)
{
  OstreamManager clout(std::cout,"setColumnTags");

  std::ofstream fout(singleton::directories().getGnuplotOutDir() + "data/" + plotFileName + ".dat",std::ios::out);

  if(fout.is_open())
  {
    unsigned int k = 0;
    for(;k < columnTags.size() -1; k++)
    {
      fout << columnTags.at(k) << _separator;
    }
    /// to prevent having a separator at the end
    fout << columnTags.at(k) << std::endl;
  }
  fout.close();
}

template< typename T >
void CSV<T>::setColumnTags(const std::vector<std::string> columnTags)
{
  setColumnTags(columnTags, _name);
}

template< typename T >
void CSV<T>::clearFile(std::string filename)
{
  /// empty the file by opening an ofstream
  std::ofstream fout(singleton::directories().getGnuplotOutDir() + "data/" + filename + ".dat",std::ios::out);
  fout.close();
}

template< typename T>
void CSV<T>::clearFile()
{
  clearFile(_name);
}


}  // namespace olb

#endif
