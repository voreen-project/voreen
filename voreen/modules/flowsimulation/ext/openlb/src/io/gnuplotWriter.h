/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Fabian Klemens
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

#ifndef GNUPLOT_WRITER_H
#define GNUPLOT_WRITER_H

#include <iomanip>
#include <iostream>
#include <vector>
#include "CSVWriter.h"
#include "CSVWriter.hh"

namespace olb {

template< typename T >
class Gnuplot {
public:

  enum AxisType {LINEAR, LOGLOG, LOGLOGINVERTED}; /// different types of data usage and axes scaling
  enum Regression {OFF, LINREG};                  /// type of Regression, off, LINear REGression, exponential Regression?,...

  /// Constructor with name for output files
  /// boolean true for real-time plotting //WARNING: experimental!
  /// Every set of paramters has its own constructor which is delegating the "large constructor"
  Gnuplot(std::string name, bool liveplot, std::string preCommand, AxisType axisType, Regression regressionType);
  explicit Gnuplot(std::string name);
  Gnuplot(std::string name, bool liveplot);
  Gnuplot(std::string name, AxisType axisType);
  Gnuplot(std::string name, AxisType axisType, Regression regressionType);
  
  /// initialises the data file
  void init();

  /// sets the data and plot file for two doubles (x and y)
  /// the plotType indicates whether the user want to plot a line graph (default = 'l')
  /// or e.g. a scatter plot ('p')
  void setData(T xValue, T yValue, std::string name = "", std::string key = "", char plotType = 'l');
  /// if no x value is given, it is just an increasing integer
  void setData(bool noXvalue, T yValue, std::string name = "", std::string key = "", char plotType = 'l');

  /// sets the data and plot file for a double and a vector of doubles (x and {y1,y2,...})
  /// for each entry of the y-axis-list, a plot Type can be specified (default: line graph {'l', 'l'})
  void setData(T xValue, std::vector<T> yValues, std::vector<std::string> names = {""}, std::string key = "right", std::vector<char> plotType = {'l','l'});
  /// if no x value is given, it is just an increasing integer
  void setData(bool noXvalue, std::vector<T> yValues, std::vector<std::string> names = {""}, std::string key = "right", const std::vector<char> plotType = {'l','l'});

  /// set labels of the plot: xLabel and yLabel
  void setLabel(std::string xLabel = "", std::string yLabel = "");

  /// writes an PDF
  void writePDF(std::string plotName = "");

  /// writes PNGs
  /// usage: first argument: numbering of png file (optional),
  /// second argument: range for the x axis (optional)
  /// third argument: name the plot in order to specify a plotID in case the
  /// user want to create more than one plot with the simulation results
  /// no arguments: writes in one file with adaptive xrange and no specific name
  void writePNG(int iT = -1, double xRange = -1, std::string plotName = "");


private:
  std::string _xLabel;
  std::string _yLabel;
  std::string _name;
  bool _liveplot;
  std::string _type;
  std::string _dataFile;
  std::string _dir;
  std::vector<std::string> _names;
  std::string _key;
  std::string _preCommand;
  AxisType _axisType; 
  Regression _regressionType; 

  /// the creation of the datafile and writing of the data in the CSV format datafile
  /// is now done by the CSV Writer. The functions setData call the CSV Writer writeDataFile function
  CSV<T> csvWriter;

  /// plotType vector stores the type of the plot the user wants to create
  /// (e.g. 'l' line graph (default) or 'p' scatterplot)
  std::vector<char> _plotTypes;
  bool _init = true;
  unsigned int _dataSize = 0;
  int _iT = -1;
  double _xRange = -1;
  T _time = 0.;

  static constexpr int _rank {0};

  /// writes a plot file for type {"plot", "png", "pdf")
  void writePlotFile(std::string type, std::string plotName = "");

protected:
  /// system command to start gnuplot (LINUX ONLY!)
  void startGnuplot(std::string plotFile, std::string plotName = "");

  /// creates the lin regression to the data
  void linRegression(std::ofstream& fout, std::string x_axisType, std::string y_axisType);

  /// scales the axes if needed
  void scaleAxes(std::ofstream& fout);
};

}  // namespace olb

#endif
