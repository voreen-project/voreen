/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013-2016 Thomas Henn, Mathias J. Krause
 *                2016-2017 Albert Mink, Maximilian Gaedtke, Markus Morhard, Mathias J. Krause
 *                2021      Nicolas Hafen, Mathias J. Krause
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

#ifndef VTK_WRITER_H
#define VTK_WRITER_H

/** This class is intended to provide the ability to write particle data.
 *  The overall structure is, however, similar to the one of the SuperVTMwriter3D
 *  rather than to the SuperParticleSysVtuWriter. This was done to merge it to the
 *  SuperVTMwriter3D someday.
 *  The following differences do not alllow this a the moment though:
 *  - The SuperVTMwriter3D is still dimension sensitiv, the present is not.
 *  - The SuperVTMwriter3D assumes a parallelized framework, the present assumes
 *    a non-parallelized framework
 *  - The Functor type differs. Those can possibly be bases on a common base class
 *    or a template argument could be provided
 *  - The SuperVTMwriter3D does not take the DESCRIPTOR as argument. The present one does
 *    as this is necessary for the underlying datatype. In order to destinguish between 
 *    particles and the lattice, the descriptor should be included though.
 *  As soon as the differences are removed, a generic 2D/3D parallized/nonparallized 
 *  VTKwriter class can be implemented providing all functionality.
 *  WARNING: Due to this, some functionality was not bothered to be implemented
 *  ( e.g. compression and binary writeout).
 *
 */

namespace olb {

//Define vtkType enumarator
enum vtkType { VTI, VTU, VTP };

//VTKwriter class
template<typename T, typename FUNCTOR, vtkType VTKTYPE>
class VTKwriter {
public:
  static constexpr bool parallel = FUNCTOR::isSuper;
  static constexpr unsigned D = FUNCTOR::d;
  /// constructor
  VTKwriter( const std::string & name, bool binary=true, bool compress=true );
  ///  put functor to _pointerVec  
  ///  to simplify writing process of several functors
  void addFunctor(FUNCTOR& f);
  ///  put functor with specific name to _pointerVec 
  ///  to simplify writing process of several functors
  void addFunctor(FUNCTOR& f, const std::string& functorName);
  ///  to clear stored functors, not yet used due to lack of necessity
  //
 



  template<bool vtmAsMaster,typename F>
  void write(int iT, std::string nameCollection, F applyFunctors);
  void write(int iT = 0);
  void write(FUNCTOR& f, int iT = 0);
  void write(std::shared_ptr<FUNCTOR> ptr_f, int iT = 0);
  ///  have to be called before calling write(int iT=0), since it creates
  //   the master pvd file, where all data files are linked!
  void createMasterFile();

protected:
  ///  performes <VTKFile ...> and <Collection>
  void preamblePVD(const std::string& fullNamePVD);
  ///  performes <DataSet timestep= ... file=namePiece />
  void dataPVD(int iT, const std::string& fullNamePVD,
                     const std::string& namePiece);
  ///  performes </Collection> and </VTKFile>
  void closePVD(const std::string& fullNamePVD);

  ///  performes <VTKFile ...> and <Collection>
  void preambleVTM(const std::string& fullNameVTM);
  /// performes </Collection> and </VTKFile>
  void closeVTM(const std::string& fullNameVTM);
  ///  performes <DataSet timestep= ... file=namePiece />
  ///  used for linking vti into pvd files
  void dataVTM(int iC, const std::string& fullNameVTM, const std::string& namePiece);
  ///  wrapper for VTM file creation
  template<bool vtmAsMaster=false>
  void writeVTM(int iT, int rankSize, std::string fileExtension, std::string nameCollection);
  ///  performes <VTKFile ...>, <ImageData ...> and <PieceExtent ...>
  void preambleVTU(const std::string& fullName, Vector<int,1> extent1);
  ///  performes </ImageData> and </VTKFile>
  void closeVTU(const std::string& fullNamePiece);
  ///  TODO: add description: connectivity, offsete, type of unscructured nodes
  void cellDataVTU(const std::string& fullName, Vector<int,1> extent1);

  ///  performes <VTKFile ...>, <ImageData ...>, <PieceExtent ...> and <PointData ...>
  void preambleVTI(const std::string& fullName,
    const LatticeR<D> extent0, const LatticeR<D> extent1,
    PhysR<T,D> origin, T delta);
  ///  performes </ImageData> and </VTKFile>
  void closeVTI(const std::string& fullNamePiece);
  ///  performes </PointData> and </Piece>
  void closePiece(const std::string& fullNamePiece);


  ///  writes points necessary for VTU
  template<unsigned sourceDim>
  void dataArrayPoints(const std::string& fullName,
                       Vector<int,sourceDim> extent1, int iC=0 );
  ///  writes functors stored at pointerVec
  template<unsigned sourceDim>
  void dataArraySingleFunctor(const std::string& fullName,
                 FUNCTOR& f,
                 Vector<int,sourceDim> extent1, int iC=0 );
private:
  mutable OstreamManager clout;
  ///  default is false, call createMasterFile() and it will be true
  bool _createFile;
  ///  determines the name of .vti and .pvd per iT
  std::string const _name;
  ///  holds added functor, to simplify the use of write function
  std::vector<FUNCTOR*> _pointerVec;
  ///  writing data base64 encoded
  bool _binary;
  ///  writing data zLib compressed
  bool _compress;
};


//VTU: vtk writer for unstructured data
template<typename T, typename DESCRIPTOR, bool parallel=true>
using VTUwriter = VTKwriter<
  T,
  typename std::conditional_t<
    parallel,
    SuperContainerF<T,DESCRIPTOR,DynamicFieldGroupsD<T,typename DESCRIPTOR::fieldList>,T>,
    ContainerF<T,DESCRIPTOR,DynamicFieldGroupsD<T,typename DESCRIPTOR::fieldList>,T>
  >,
  VTU
>;

//VTI (image) version, corresponding to original SuperVTMwriterXD
template<typename T, typename DESCRIPTOR>
using VTIwriter = VTKwriter<
  T,
  SuperF3D<T,T>,
  VTI
>;

}  // namespace OLB

#endif
