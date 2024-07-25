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

#ifndef GNUPLOT_WRITER_HH
#define GNUPLOT_WRITER_HH

#include <fstream>
#include <iostream>
#include <unistd.h>
#include "core/singleton.h"
#include "io/fileName.h"
#include "gnuplotWriter.h"
#include "utilities/vectorHelpers.h"

namespace olb {
/// Constructor with name of outputFiles
/// boolean true for real-time plotting //WARNING: experimental!
template< typename T >
Gnuplot<T>::Gnuplot(std::string name, bool liveplot, std::string preCommand, AxisType axisType, Regression regressionType) 
  : _name(name),
    _liveplot(liveplot),
    _dataFile(singleton::directories().getGnuplotOutDir()+"data/"+_name+".dat"),
    _dir(singleton::directories().getGnuplotOutDir()),
    _preCommand(preCommand),
    _axisType(axisType),
    _regressionType(regressionType),
    csvWriter(name)
{ }

/// overload the constructor using delegating constructors. 
/// This has the advantage that there is no necessity for default values which might cause some trouble
template <typename T>
Gnuplot<T>::Gnuplot(std::string name) : Gnuplot(name, false, "", LINEAR, OFF) {}

template< typename T >
Gnuplot<T>::Gnuplot(std::string name, bool liveplot) : Gnuplot(name, liveplot,"", LINEAR, OFF) {}

template <typename T>
Gnuplot<T>::Gnuplot(std::string name, AxisType axisType) : Gnuplot(name, false, "", axisType, OFF) {}

template< typename T >
Gnuplot<T>::Gnuplot(std::string name, AxisType axisType, Regression regressionType) : Gnuplot(name, false,"", axisType, regressionType) {}


/// writes the data and plot file for two doubles (x and y)
/// plotType indicates whether you want a linegraph 'l' (default) or a scatterplot 'p' (default: 'l')
template< typename T >
void Gnuplot<T>::setData(T xValue, T yValue, std::string name, std::string key, char plotType)
{
  if (_init) {
    _dataSize = 1;
    _key = key;
    _names = {name};
    _plotTypes = {plotType};

    if (_liveplot) {
      writePlotFile("plot");
    }
  }
  csvWriter.writeDataFile(xValue, yValue);

  if (_liveplot && _init) {
    startGnuplot("plot");
  }

  _init = false;
  return;
}

/// writes the data and plot file for two doubles (x and y), where x is increasing integer
template< typename T >
void Gnuplot<T>::setData(bool noXvalue, T yValue, std::string name, std::string key, char plotType)
{
  T xValue = _time;
  setData(xValue, yValue, name, key, {plotType});
  _time++;
}

/// writes the data and plot file for a double and a vector of doubles (x and y1,y2,...)
/// plotType indicates whether you want a linegraph 'l' (default) or a scatterplot 'p': (default: {'l','l'})
/// The position in the vector 'plotType'{'l', 'p'} is linked to the rank of the y-axis (y1, y2) :
/// y1 is plotted in form of a line plot & y2 is plotted in form of a scatterplot
template< typename T >
void Gnuplot<T>::setData(T xValue, std::vector<T> yValues, std::vector<std::string> names, std::string key, std::vector<char> plotType)
{
  if (_init) {
    _dataSize = yValues.size();
    _key = key;
    _names = names;
    _plotTypes = plotType;
    if (_names.size() < _dataSize) {
      for (unsigned int i = _names.size(); i < _dataSize; i++) {
        _names.push_back("");
      }
    }
    if (_plotTypes.size() < _dataSize) {
      for (unsigned int i = _plotTypes.size(); i < _dataSize; i++) {
        _plotTypes.push_back('l');
      }
    }
    if (_liveplot) {
      writePlotFile("plot");
    }
  }
  csvWriter.writeDataFile(xValue,yValues);

  if (_liveplot && _init) {
    startGnuplot("plot");
  }

  _init = false;
  return;
}

/// create lin Regression to the given Datasets in the given scaling, e.g. loglog, logloginverted, linear,....
/// The kind of scaling is provided by the parameters xAxisType and yAxisType
template<typename T>
void Gnuplot<T>::linRegression(std::ofstream& fout, std::string xAxisType, std::string yAxisType)
{ 
  fout << "set fit quiet\n";

  for (unsigned int i = 0; i < _dataSize; ++i) {
    fout << "f" << i+1 << "(x) = m" << i+1 << " * x + b"<< i+1 << "\n"; //Match f2 to the dataset 1:3 with f2(x) = m2 * x + b2 and so on
    fout << "fit f" << i+1 << "(x) path using ("<< xAxisType << "$1)):(" << yAxisType << "$" << i+2 << ")) via m" << i+1 << ", b" << i+1 <<"\n";
  }
  fout << "\n";

  ///plotting the data and for each set 1:2 / 1:3,... the lin regression is plotted
  ///IMPORTANT the word "plot" must be included in every kind of regression used in the context auf the plotting
  fout << "plot ";
  for (unsigned int i = 0; i < _dataSize; ++i) {
    fout << "f"<< i+1 <<"(x) title sprintf('regr " << _names[i] << ", gradient = %.3f', m" << i+1 << ") lc " << i+1 << " axis x1y1,";
  }

  return;
}

/// scales the axes if necessary
/// to add new kinds of scaling, just add another case, the labels for the x and y axes will be expanded by the kind of scaling (if it's other than linear)
/// this labeling provides a labeling of the axes even if the function setLabel isn't specifically called 
template<typename T>
void Gnuplot<T>::scaleAxes(std::ofstream& fout)
{
  switch (_axisType)
  {
    case (LOGLOG): ///Difference between LOGLOG and LOGLOGINVERTED is just 1/N in the x-axes
    case(LOGLOGINVERTED):
    {
      fout << "set xtics nomirror" << "\n";
      fout << "set ytics nomirror" << "\n";

      fout << "set format y '10^{%.2f}'" << "\n"; ///this command will change the tics-labels, the {%.2f} is replaced by the original tic-label using two
      fout << "set format x '10^{%.2f}'" << "\n"; ///digits after the comma
    } break;

    case LINEAR:
    default:
    {} break;
  }
  
  fout << "\n";
  return;
}


/// writes the data and plot file for a double and a vector of doubles (x and y1,y2,...), where x is increasing integer
template< typename T >
void Gnuplot<T>::setData(bool noXvalue, std::vector<T> yValues, std::vector<std::string> names, std::string key, std::vector<char> plotType)
{
  T xValue = _time;
  setData(xValue, yValues, names, key, plotType);
  _time++;
}


/// writes an PDF
template< typename T >
void Gnuplot<T>::writePDF(std::string plotName)
{
  if (!_init) {
    writePlotFile("pdf", plotName);
    startGnuplot("plotPDF", plotName);
  }
  return;
}


/// writes PNGs
/// usage: first argument: numbering of png file
/// second argument: range for the x axis
/// thrid argument: specifies the name of the plot in case the user wants to
/// create more than one plot with the simulation results (default: plotName = "")
/// no arguments: writes consecutive numbers with adaptive xrange
template< typename T >
void Gnuplot<T>::writePNG(int iT, double xRange, std::string plotName)
{
  if (!_init) {
    _iT = iT;
    _xRange = xRange;

    /// initialize the writePlotFile for Gnuplot with the type and the name of the output data
    writePlotFile("png", plotName);
    startGnuplot("plotPNG", plotName);
  }
  return;
}

/// This function is kind of the "heart" of the data analysis as it writes the gnuplot file that will be executed
/// plotName specifies the name of the plot in case the user wants to create more than
/// one plot with the simulation results (default: plotName = "")
template< typename T >
void Gnuplot<T>::writePlotFile(std::string type, std::string plotName)
{
  if (singleton::mpi().getRank() == _rank ) {
    std::ofstream fout;

    std::string plotFile;
    if (_liveplot && type == "plot") {
      plotFile = singleton::directories().getGnuplotOutDir()+"data/plot.p";
    }
    else if (type == "pdf") {
      plotFile = singleton::directories().getGnuplotOutDir()+"data/plotPDF"+plotName+".p";
    }
    else if (type == "png") {
      plotFile = singleton::directories().getGnuplotOutDir()+"data/plotPNG"+plotName+".p";
    }
    else {
      std::cout << "WARNING: invalid Gnuplot type={'', 'plot'; 'pdf', 'png'}" << std::endl;
      exit(-1);
    }

    fout.open(plotFile.c_str(), std::ios::trunc);
    fout << "set key " << _key << "\n";

    if (type=="pdf") {
      fout << "set terminal pdf enhanced" << "\n"
           << "set output '"<<_dir<<_name<<".pdf'" << "\n";
    }
    if (type=="png") {
      if ( !util::nearZero(_xRange+1) ) {
        fout << "set xr[0:"<< _xRange <<"]" << "\n";
      }
      fout << "set terminal png" << "\n"
           << "set output '"<<_dir<<_name;
      if (_iT != -1) {
        fout <<"_"<<_iT;
      }
      fout <<".png'" << "\n";
    }
    
    /// set precommands
    fout << _preCommand << "\n";

    
    /// Scaling the axes if necessary
    scaleAxes(fout); 
    
    /// Labelling the axes

    fout << "set xlabel '" << _xLabel << "'" << "\n";
    fout << "set ylabel '" << _yLabel << "'" << "\n";



    /// These variables will be placed directly in the plot-command
    std::string xAxisType; /// Variable to decide how the data should be used and is scaled (e.g. log10(data), ...)
    std::string yAxisType; /// Differentiation between x and y axis is necessary to implement e.g. log10(1/$1) instead of log10($1)

    switch (_axisType)
    {
      case LOGLOG:
      {
        xAxisType = "log10(";
        yAxisType = "log10(";
      } break;

      case LOGLOGINVERTED:
      {
        xAxisType = "log10(1/";
        yAxisType = "log10(";
      } break;

      case LINEAR:
      default:
      {
        xAxisType = "(";
        yAxisType = "(";
      }
    }

    fout << "path = '" << _dir << "data/"<<_name << ".dat'" << "\n"; /// Path to the data, saves a lot of code in the gnuplot file and within this file

    

    /// The necessary statements for a regression will be added if needed
    /// IMPORTANT for later additions: the word "plot" needs to be included in every case separatly, see e.g. linRegression
    /// One must also add the buzzword (like LINREG, OFF) in the enum declaration in gnuplotWriter.h

    /// there might occur the error notification "not data to fit" while the compiler is running which doesn't directly indicate whether
    /// the plot/regression was succesfull or not
    switch (_regressionType)
    {
      case LINREG:
      {  /// LINear REGression of the datasets
        linRegression(fout, xAxisType, yAxisType);  /// Create the lin Regression of the datasets, includes the plotting of the regression AND the data
      } break;
      
      case OFF: /// No Regression, just plot the data, OFF same as default
      default:
      {
        fout << "plot ";
      } break;

    }
    /// vector which holds the information about the plotType
    /// (e.g. scatterplot 'p' or lineplot 'l': default {'l','l'}) 

    /// plotting the data and for each set 1:2 / 1:3,... 

    unsigned int i = 0;
    for ( ; i < _dataSize - 1; ++i) {
      fout << "path u (" << xAxisType << "$1)):(" << yAxisType << "$" << i+2 << ")) w " << _plotTypes[i] << " t '" << _names[i] << "' lc " << i+1 << " axis x1y1,";
    }
    fout << "path u (" << xAxisType << "$1)):(" << yAxisType << "$" << i+2 << ")) w " << _plotTypes[i] << " t '" << _names[i] << "' lc " << i+1 << " axis x1y1";

    fout << "\n";
    if (_liveplot && type=="plot") {
      fout << "pause -1" << "\n"
          << "reread" << "\n";
    }
    fout.close();

    
  }
  return;
}


/// set Label of the gnuplotPlot; xLabel and yLabel
template< typename T >
void Gnuplot<T>::setLabel(std::string xLabel, std::string yLabel)
{
  _xLabel = xLabel;
  _yLabel = yLabel;
}


/// system command to start gnuplot (LINUX ONLY!)
/// plotName indicates the name of the plot in case the user wants to create
/// more than one plot with the simulation results (default: plotName = "")
template< typename T >
void Gnuplot<T>::startGnuplot(std::string plotFile, std::string plotName)
{
#ifdef WIN32
  std::cout << "GNUPLOT WORKS ONLY WITH LINUX" << std::endl;
//  exit (EXIT_FAILURE);
  return;
#endif
#ifndef WIN32
  if (singleton::mpi().getRank() == _rank) {
    if (!system(nullptr)) {
      exit (EXIT_FAILURE);
    }
    const std::string command = "gnuplot -persistent "+_dir+"data/"+plotFile+plotName+".p > /dev/null &";
    if ( system(command.c_str()) ) {
      std::cout << "Error at GnuplotWriter" << std::endl;
    }
  }
  return;
#endif
}
}  // namespace olb

#endif
