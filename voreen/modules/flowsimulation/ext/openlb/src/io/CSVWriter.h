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

#ifndef CSV_WRITER_H
#define CSV_WRITER_H

#include <iomanip>
#include <iostream>
#include <vector>

namespace olb {

template< typename T >
class CSV {
public:

  /// Constructor with name for output file
  CSV(std::string name, char separator, std::vector<std::string> columnTags);
  explicit CSV(std::string name = "unnamed");
  CSV(std::string name, char separator);
  CSV(std::string name, std::vector<std::string> columnTags);

  /// former datFileOut functions
  /// these functions create a datafile in the csv format
  /// in order to write in a separated datafile as the one named in the constructor,
  /// just put the name as the plotFileName parameter
  void writeDataFile(T xValue, T yValue, const std::string& plotFileName, int precision = 16);
  void writeDataFile(T xValue, const std::vector<T>& yValues, const std::string& plotFileName, int precision = 16);
  /// writes the data file for two doubles (x and y)
  void writeDataFile(T xValue, T yValue, int precision = 16);
  /// writes the data file for one double and a vector of doubles (x and y1,y2,...)
  void writeDataFile(T xValue,const std::vector<T>& yValues, int precision = 16);

  /// adds column tags at the beginning of the csv data file.
  /// saves the content of the file and rewrites it with the tags at the top and then inserting the content again
  void setColumnTags(const std::vector<std::string> columnTags, std::string& plotFileName);
  void setColumnTags(const std::vector<std::string> columnTags);

  /// clears the file
  void clearFile(std::string filename);
  void clearFile();

private:
  std::string _name;
  std::string _dataFile;
  std::string _dir;
  char _separator;
  std::vector<std::string> _columnTags;

  static constexpr int _rank {0};
};


}  // namespace olb

#endif
