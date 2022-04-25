/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Julius Jessberger, Adrian Kummerlaender
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


#ifndef SOLVER_PARAMETERS_H
#define SOLVER_PARAMETERS_H


#include <iostream>

#include "io/ostreamManager.h"
#include "io/xmlReader.h"
#include "utilities/benchmarkUtil.h"
#include "names.h"
#include "utilities/typeMap.h"

namespace olb {

namespace parameters {

struct ParameterBase {
  virtual void initialize() { }
};

// ------------------------------ Simulation ----------------------------------

/** Base struct to keep the parameters that are necessary for the simulation.
 * Class supports reading simulation parameters from input files (e.g. XML) as
 * well as an interface for other parts of the program.
 * This class is intended to be derived from: at least, a UnitConverter has to
 * be provided additionally.
 */
template<typename T>
struct SimulationBase : public ParameterBase
{
  using BT = BaseType<T>;

  BT                                    startUpTime {0};
  BT                                    maxTime {1};
  BT                                    physBoundaryValueUpdateTime {0.1};

  bool                                  pressureFilter {false};

  // no of cuboids per mpi process
  int                                   noC     {1};
  int                                   overlap {3};

  virtual void initialize() { }
};

/// All the simulation parameters are read directly from an xml file
template<typename T, typename LATTICES>
struct XmlSimulation : public SimulationBase<T> {

  XMLreader     const* xml;

  using NAME = typename LATTICES::keys_t::template get<0>;
  using descriptor = typename LATTICES::template value<NAME>;
  std::shared_ptr<UnitConverter<T,descriptor> const> converter;
};


// ------------------------------ Stationarity --------------------------------

struct StationarityBase : public ParameterBase
{ };

/** All parameters that are necessary for checking whether the simulation
 * became stationary.
 * Template parameters are the arithmetic type T and (optionally) the names of
 * lattices for which stationarity is to be checked. If none are given,
 * NavierStokes is set as a default.
 */
template<typename T, typename... STAT_LATTICES>
struct Stationarity : public StationarityBase
{
  constexpr static unsigned numberOfArguments = sizeof...(STAT_LATTICES);
  // list of lattices for which stationarity is to be checked
  using stat_lattices = std::conditional_t <(numberOfArguments == 0),
    meta::list<names::NavierStokes>,
    meta::list<STAT_LATTICES...>>;
  constexpr static unsigned numberOfStationaryLattices = stat_lattices::size;

  enum ConvergenceType {MaxLatticeVelocity, AverageEnergy, AverageRho};

  std::array <ConvergenceType,numberOfStationaryLattices> convergenceType {MaxLatticeVelocity};
  std::array <T,numberOfStationaryLattices>               physInterval {1};
  std::array <T,numberOfStationaryLattices>               epsilon {1e-3};

  Stationarity() = default;

  Stationarity(ConvergenceType type, T interval_, T epsilon_)
  {
    std::fill(convergenceType.begin(), convergenceType.end(), type);
    std::fill(physInterval.begin(), physInterval.end(), interval_);
    std::fill(epsilon.begin(), epsilon.end(), epsilon_);
  }

  Stationarity(
    std::array <ConvergenceType,numberOfStationaryLattices> type,
    std::array <T,numberOfStationaryLattices> interval_,
    std::array <T,numberOfStationaryLattices> epsilon_)
   : convergenceType(type), physInterval(interval_), epsilon(epsilon_)
  { }
};


// ------------------------------ Output --------------------------------------

/** Struct to keep parameters which characterize the output.
 * So far, bash/ gnuplot/ image/ vtk are supported.
 */
template<typename T, typename LatticeLog= names::NavierStokes>
struct OutputBase : public ParameterBase
{
  using BT = BaseType<T>;
  std::string                              name         {"unnamed"};
  std::string                              olbDir       {"../../../"};
  std::string                              outputDir    {"./tmp/"};

  bool                                     verbose            {true};
  bool                                     printLogConverter  {true};

  BT                                       logT;
  using printLatticeStatistics =           meta::list<LatticeLog>;

  int                                      timerPrintMode     {2};

  bool                                     gnuplotOutput      {false};
  std::string                              gnuplotFilename    {"unnamed"};
  BT                                       gnuplotSaveT;
  bool                                     imageOutput        {false};
  std::string                              imageFilename      {"unnamed"};
  BT                                       imageSaveT;
  bool                                     vtkOutput          {false};
  std::string                              vtkFilename        {"unnamed"};
  BT                                       vtkSaveT;

  OutputBase() = default;

  OutputBase(
    std::string name_,
    std::string olbDir_,
    std::string outputDir_,
    bool verbose_,
    bool printLogConverter_,
    BT logT_,
    int timerPrintMode_,
    bool gnuplotOutput_,
    std::string gnuplotFilename_,
    BT gnuplotSaveT_,
    bool imageOutput_,
    std::string imageFilename_,
    BT imageSaveT_,
    bool vtkOutput_,
    std::string vtkFilename_,
    BT vtkSaveT_
  ) : name(name_), olbDir(olbDir_), outputDir(outputDir_), verbose(verbose_),
    printLogConverter(printLogConverter_), logT(logT_),
    timerPrintMode(timerPrintMode_), gnuplotOutput(gnuplotOutput_),
    gnuplotFilename(gnuplotFilename_), gnuplotSaveT(gnuplotSaveT_),
    imageOutput(imageOutput_), imageFilename(imageFilename_),
    imageSaveT(imageSaveT_), vtkOutput(vtkOutput_), vtkFilename(vtkFilename_),
    vtkSaveT(vtkSaveT_)
  { }

  virtual void initialize() override
  {
    singleton::directories().setOlbDir(olbDir);
    singleton::directories().setOutputDir(outputDir);
  }
};


// ------------------------------ Results -------------------------------------

/** Struct to keep results of the simulation in order to provide communication
 * with other parts of the program.
 * This class is intended to be derived from.
 */
struct ResultsBase : public ParameterBase { };


// ------------------ Parameter reading functionality ------------------------

template<typename PARAMETERS>
struct ReaderBase
{
  mutable OstreamManager clout {std::cout, "ParameterReader"};
  std::shared_ptr<PARAMETERS> params;

  ReaderBase(std::shared_ptr<PARAMETERS> params_) : params(params_)
  { }
};

/** Base struct for reading parameters from files.
 * So far, reading from xml is supported. If you want to let your own parameter
 * struct read from xml, you may realize this as template specializations,
 * similar to the syntax in the examples below. The default (with no
 * specialization) does nothing.
 * If you want to add your own interface, you may add read(...)-methods to the
 * given structures which read the parameters from your file type.
 */
template<typename PARAMETERS>
struct Reader : public ReaderBase<PARAMETERS> {
  Reader(std::shared_ptr<PARAMETERS> params_) : Reader::ReaderBase(params_)
  { }

  virtual void read(XMLreader const& xml) { }
};

template<typename PARAMETERS>
Reader(std::shared_ptr<PARAMETERS> params_) -> Reader<PARAMETERS>;


template<typename T>
struct Reader<OutputBase<T>> : public ReaderBase<OutputBase<T>>
{
  using ReaderBase<OutputBase<T>>::ReaderBase;

  void read(XMLreader const& xml)
  {
    // Filenames
    xml.readOrWarn<std::string>(this->clout, "Application", "Name", "", this->params->name, true, true, true);
    xml.readOrWarn<std::string>(this->clout, "Application", "OlbDir", "", this->params->olbDir, true, true, true);
    xml.readOrWarn<std::string>(this->clout, "Output", "OutputDir", "", this->params->outputDir, true, true, true);

    // For Console Output
    xml.readOrWarn<bool>(this->clout, "Output", "Log", "VerboseLog", this->params->verbose, true, false, true);
    xml.readOrWarn<bool>(this->clout, "Output", "PrintLogConverter", "", this->params->printLogConverter, true, false, false);
    xml.readOrWarn<BaseType<T>>(this->clout, "Output", "Log", "SaveTime", this->params->logT, true, false, true);
    xml.readOrWarn<int>(this->clout, "Output", "Timer", "PrintMode", this->params->timerPrintMode, true, false, false);

    // For Gnuplot Output
    xml.readOrWarn<bool>(this->clout, "Output", "VisualizationGnuplot", "Output", this->params->gnuplotOutput, true, false, true);
    if (this->params->gnuplotOutput) {
      this->params->gnuplotFilename = this->params->name;
      xml.readOrWarn<std::string>(this->clout, "Output", "VisualizationGnuplot", "Filename", this->params->gnuplotFilename, true, false, false);
      xml.readOrWarn<BaseType<T>>(this->clout, "Output", "VisualizationGnuplot", "SaveTime", this->params->gnuplotSaveT, true, false, true);
    }

    // For Image Output
    xml.readOrWarn<bool>(this->clout, "Output", "VisualizationImages", "Output", this->params->imageOutput, true, false, true);
    if (this->params->imageOutput) {
      this->params->imageFilename = this->params->name;
      xml.readOrWarn<std::string>(this->clout, "Output", "VisualizationImages", "Filename", this->params->imageFilename, true, false, false);
      xml.readOrWarn<BaseType<T>>(this->clout, "Output", "VisualizationImages", "SaveTime", this->params->imageSaveT, true, false, true);
    }

    // For VTK Output
    xml.readOrWarn<bool>(this->clout, "Output", "VisualizationVTK", "Output", this->params->vtkOutput, true, false, true);
    if (this->params->vtkOutput) {
      this->params->vtkFilename = this->params->name;
      xml.readOrWarn<std::string>(this->clout, "Output", "VisualizationVTK", "Filename", this->params->vtkFilename, true, false, false);
      xml.readOrWarn<BaseType<T>>(this->clout, "Output", "VisualizationVTK", "SaveTime", this->params->vtkSaveT, true, false, true);
    }
  }
};

template<typename T>
struct Reader<SimulationBase<T>> : public ReaderBase<SimulationBase<T>>
{
  using ReaderBase<SimulationBase<T>>::ReaderBase;

  void read(XMLreader const& xml)
  {
    using BT = BaseType<T>;

    xml.readOrWarn<BT>(this->clout, "Application", "PhysParameters", "StartUpTime", this->params->startUpTime, false, false, false);
    xml.readOrWarn<BT>(this->clout, "Application", "PhysParameters", "PhysMaxTime", this->params->maxTime, false, true, true);
    xml.readOrWarn<int>(this->clout, "Application", "Mesh", "noCuboidsPerProcess", this->params->noC, true, false, true);
    xml.readOrWarn<bool>(this->clout, "Application", "PressureFilter", "", this->params->pressureFilter, true, false, false);

    // TODO: Find a better place for dependent default values.
    this->params->physBoundaryValueUpdateTime = this->params->maxTime / BT(100); // 1% of max. time as default value
    xml.readOrWarn<BT>(this->clout, "Application", "PhysParameters", "BoundaryValueUpdateTime", this->params->physBoundaryValueUpdateTime, true, false, true);
  }
};

// todo: enable reading for multiple stationary lattices
template<typename T>
struct Reader<Stationarity<T>> : public ReaderBase<Stationarity<T>>
{
  using ReaderBase<Stationarity<T>>::ReaderBase;

  void read(XMLreader const& xml)
  {
    std::string type = "MaxLatticeVelocity";
    xml.readOrWarn<std::string>(this->clout, "Application", "ConvergenceCheck", "Type", type, true, false, true);
    this->params->convergenceType[0] = (type == "AverageEnergy") ? Stationarity<T>::AverageEnergy : Stationarity<T>::MaxLatticeVelocity;
    xml.readOrWarn<BaseType<T>>(this->clout, "Application", "ConvergenceCheck", "Interval", this->params->physInterval[0], true, false, true);
    xml.readOrWarn<BaseType<T>>(this->clout, "Application", "ConvergenceCheck", "Residuum", this->params->epsilon[0], true, false, true);
  }
};

template<typename T, typename LATTICES>
struct Reader<XmlSimulation<T,LATTICES>> : public ReaderBase<XmlSimulation<T,LATTICES>>
{
  using ReaderBase<XmlSimulation<T,LATTICES>>::ReaderBase;

  void read(XMLreader const& xml_)
  {
    Reader<SimulationBase<T>>(this->params).read(xml_);
    this->params->xml = &xml_;

    using NAME = typename LATTICES::keys_t::template get<0>;
    using descriptor = typename LATTICES::template value<NAME>;
    this->params->converter = std::shared_ptr<UnitConverter<T,descriptor>>(createUnitConverter<T,descriptor>(xml_));
  }
};

} // namespace parameters


// ------------------ Creator functions for XML interface ---------------------

template <typename PARAMETERS>
std::shared_ptr<PARAMETERS> createParameters(XMLreader const& xml)
{
  auto result = std::make_shared<PARAMETERS>();
  parameters::Reader<PARAMETERS>(result).read(xml);
  result->initialize();
  return result;
}


} // namespace olb

#endif
