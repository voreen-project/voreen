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

#ifndef VTK_WRITER_HH
#define VTK_WRITER_HH



namespace olb {


template<typename T, typename FUNCTOR, vtkType VTKTYPE>
VTKwriter<T,FUNCTOR,VTKTYPE>::VTKwriter( const std::string & name, bool binary, bool compress )
  : clout( std::cout, "VTKwriter" ), _createFile(false), _name(name), _binary(binary), _compress(compress)
{}


template<typename T, typename FUNCTOR, vtkType VTKTYPE>
void VTKwriter<T,FUNCTOR,VTKTYPE>::addFunctor(FUNCTOR& f)
{
  _pointerVec.push_back(&f);
}

template<typename T, typename FUNCTOR, vtkType VTKTYPE>
void VTKwriter<T,FUNCTOR,VTKTYPE>::addFunctor(FUNCTOR& f, const std::string& functorName)
{
  f.getName() = functorName;
  _pointerVec.push_back(&f);
}


template<typename T, typename FUNCTOR, vtkType VTKTYPE>
template<bool vtmAsMaster>
void VTKwriter<T,FUNCTOR,VTKTYPE>::writeVTM( int iT, int rankSize, std::string fileExtension, std::string nameCollection)
{
  //Evaluate vtmAsMaster
  std::string nameVTM;
  std::string pathPrefixFileData;
  if constexpr(!vtmAsMaster){
    //Name of vtm
    nameVTM =  "data/" + createFileName(nameCollection, iT ) + ".vtm";
    //path prefix
    pathPrefixFileData = "";
    //Create PVD
    std::string fullNamePVD = singleton::directories().getVtkOutDir()
                            + createFileName(nameCollection) + ".pvd";
    dataPVD(iT, fullNamePVD, nameVTM);
  } else {
    nameVTM = createFileName(nameCollection, iT) + ".vtm";
    pathPrefixFileData = "data/";
  }
  //Create full name VTM
  std::string fullNameVTM = singleton::directories().getVtkOutDir() + nameVTM;
  //Write preamble VTM
  preambleVTM(fullNameVTM);
  //Loop over all cuboids (globIcs) 
  for (int iC = 0; iC < rankSize; iC++) {
    std::string nameFileData = pathPrefixFileData + createFileName(nameCollection, iT, iC) + fileExtension;
    // puts name of .vti piece to a .pvd file [fullNameVTM]
    dataVTM(iC, fullNameVTM, nameFileData);
    // adds a namePiece to master.pvd file.
  }
  // To do so we overwrite closePVD() and add new entry.
  closeVTM(fullNameVTM);            // timestep
}


template<typename T, typename FUNCTOR, vtkType VTKTYPE>
template<bool vtmAsMaster, typename F>
void VTKwriter<T,FUNCTOR,VTKTYPE>::write(int iT, std::string nameCollection, F applyFunctors)
{

  std::string fileExtension;

  int rank = 0;
  int rankSize;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif

  if ( _pointerVec.empty() ) {
    clout << "Error: Did you add a Functor ?" << std::endl;
  } else {

    //PARALLEL FUNCTORTYPE
    if constexpr (parallel) {

      //VTI
      if constexpr(VTKTYPE==VTI){
        fileExtension = ".vti";

        //Retrieve cuboidGeometry from first functor (TODO: include into common functortype)
        auto& superStructure = dynamic_cast<SuperF<D,T>*>(_pointerVec[0])->getSuperStructure();
        auto const& cGeometry = superStructure.getCuboidGeometry();
        LoadBalancer<T> const& loadBalancer = superStructure.getLoadBalancer();
        rankSize = loadBalancer.getRankSize();
        //Ensure no gaps between vti files (TODO: still necessary?)
        superStructure.communicate();
        //Retrieve delta
        const T delta = cGeometry.getMotherCuboid().getDeltaR();

        //Cor each cuboid
        for (int iCloc = 0; iCloc < loadBalancer.size(); iCloc++) { 

          //Retrieve global IC
          int globIC = loadBalancer.glob(iCloc);
          //Retrieve BlockGeometry
          auto blockGeometry = cGeometry.get(globIC);
          //Retrieve extent
          Vector<int,D> extent0(-1);
          Vector<int,D> extent1( blockGeometry.getExtent() );

          //Name of Data file
          std::string fullNameFileData = singleton::directories().getVtkOutDir()
                                    + "data/" + createFileName(nameCollection, iT, globIC ) + fileExtension;

          // get dimension/extent for each cuboid
          LatticeR<D> originCuboid(0.);
          T originPhysR[D] = {T()};
          cGeometry.getPhysR(originPhysR, globIC, originCuboid);
          //Write preamble VTU
          preambleVTI(fullNameFileData, extent0, extent1, originPhysR, delta);
          //Loop over functors or call individual one
          applyFunctors(fullNameFileData, extent1, globIC);
          //Close piece and vti
          closePiece(fullNameFileData);
          closeVTI(fullNameFileData);
        }

      //VTU
      } else if constexpr(VTKTYPE==VTU){
        fileExtension = ".vtu";

        //Retrive loadBalancer from first functor (TODO: include into common functortype)
        LoadBalancer<T> const& loadBalancer = _pointerVec[0]->getLoadBalancer();
        rankSize = loadBalancer.getRankSize();

        //Cor each cuboid
        for (int iCloc = 0; iCloc < loadBalancer.size(); iCloc++) { 
          //Retrieve global IC
          int globIC = loadBalancer.glob(iCloc);
          // Retrieve container size (number of particles) (TODO: include into common functortype)
          Vector<int,1> extent1 = _pointerVec[0]->getContainerF(iCloc).getContainerSize();
          //Name of VTU
          std::string fullNameFileData = singleton::directories().getVtkOutDir()
                                    + "data/" + createFileName(nameCollection, iT, globIC ) + fileExtension;
          //Write preamble VTU
          preambleVTU(fullNameFileData,extent1);
          //Loop over functors or call individual one
          applyFunctors(fullNameFileData, extent1, globIC);
          // Write celldata (connectivity, offset, types)
          this->cellDataVTU(fullNameFileData, extent1);
          // Write points 
          this->dataArrayPoints(fullNameFileData, extent1, globIC);
          // CloseVTU
          closeVTU(fullNameFileData);
        }

      //VTP
      } else {
        fileExtension = ".vtp";
        clout << "Error: VTP type not implemented yet" << std::endl;
      }

      //Write VTM (only master)
      if (rank == 0) { writeVTM<vtmAsMaster>(iT,rankSize,fileExtension,nameCollection); }

    //UNPARRALLIZED FUNCTORTYPE (for now only relevant for VTU)
    } else {
      fileExtension = ".vtu";

      //if master
      if (rank == 0) {
        // Full name pvd
        std::string fullNamePVD = singleton::directories().getVtkOutDir()
                                  + createFileName(nameCollection) + ".pvd";
        // Retrieve container size (number of particles)(TODO: include into common functortype)
        Vector<int,1> extent1(_pointerVec[0]->getContainerSize());
        // Name VTU
        std::string nameFileData =  "data/" + createFileName(nameCollection, iT ) + fileExtension;
        // Create PVD 
        dataPVD(iT, fullNamePVD, nameFileData);
        // Full Name VTU
        std::string fullNameFileData = singleton::directories().getVtkOutDir() + nameFileData;
        // Write preample VTU
        preambleVTU(fullNameFileData, extent1);
        //Loop over functors or call individual one
        applyFunctors(fullNameFileData, extent1);
        // Write celldata (connectivity, offset, types)
        this->cellDataVTU(fullNameFileData, extent1);
        // Write points 
        this->dataArrayPoints(fullNameFileData, extent1);
        // CloseVTU
        closeVTU(fullNameFileData);
      } // master only
    } // if parralized type
  } //if ( _pointerVec.empty() )
}


template<typename T, typename FUNCTOR, vtkType VTKTYPE>
void VTKwriter<T,FUNCTOR,VTKTYPE>::write(int iT)
{
  //VTI
  if constexpr(VTKTYPE==VTI){
    //Define functors to be applied
    auto applyFunctors = [&](std::string fullNameFileData, Vector<int,D> extent1, int globIC=0 )
    {
      for (int iF=0; iF<_pointerVec.size(); ++iF) {
        this->dataArraySingleFunctor(fullNameFileData, *_pointerVec[iF], extent1, globIC);
      }
    };
    //Call write
    write<false>(iT,_name,applyFunctors);

  //VTU
  } else if constexpr(VTKTYPE==VTU){
    //Define functors to be applied
    auto applyFunctors = [&](std::string fullNameFileData, Vector<int,1> extent1, int globIC=0 )
    {
      for (std::size_t iF=1; iF<_pointerVec.size(); ++iF) {
        this->dataArraySingleFunctor(fullNameFileData, *_pointerVec[iF], extent1, globIC);
      }
    };
    //Call write
    write<false>(iT,_name, applyFunctors);
  }
}


template<typename T, typename FUNCTOR, vtkType VTKTYPE>
void VTKwriter<T,FUNCTOR,VTKTYPE>::write(FUNCTOR& f, int iT)
{ 
  //VTI
  if constexpr(VTKTYPE==VTI){
    //Define functor
    auto applyFunctors = [&](std::string fullNameFileData, Vector<int,D> extent1, int globIC=0 ){
      this->dataArraySingleFunctor(fullNameFileData, f, extent1, globIC); };
    //Call write
    write<true>(iT,f.getName(),applyFunctors);

  //VTU
  } else if constexpr(VTKTYPE==VTU){
    //Define functor
    auto applyFunctors = [&](std::string fullNameFileData, Vector<int,1> extent1, int globIC=0 ){
      this->dataArraySingleFunctor(fullNameFileData, f, extent1, globIC); };
    //Call write
    write<true>(iT,f.getName(),applyFunctors);
  }

}

template<typename T, typename FUNCTOR, vtkType VTKTYPE>
void VTKwriter<T,FUNCTOR,VTKTYPE>::write(std::shared_ptr<FUNCTOR> ptr_f, int iT)
{ 
  write(*ptr_f, iT);
}


template<typename T, typename FUNCTOR, vtkType VTKTYPE>
void VTKwriter<T,FUNCTOR,VTKTYPE>::createMasterFile()
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif
  if (rank == 0) {
    std::string fullNamePVD = singleton::directories().getVtkOutDir()
                                    + createFileName(_name) + ".pvd";
    preamblePVD(fullNamePVD);
    closePVD(fullNamePVD);
    _createFile = true;
  }
}

template<typename T, typename FUNCTOR, vtkType VTKTYPE>
void VTKwriter<T,FUNCTOR,VTKTYPE>::preambleVTU(
  const std::string& fullName, Vector<int,1> extent1)
{
  std::ofstream fout(fullName, std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }
  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" ";
  if (_compress) {
    fout << "byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
  }
  else {
    fout << "byte_order=\"LittleEndian\">\n";
  }
  fout << "<UnstructuredGrid>" << std::endl;
  fout << "<Piece NumberOfPoints=\"" << extent1[0]
       << "\" NumberOfCells=\"" << extent1[0] << "\">"
       << std::endl;
  fout << "<PointData Vectors=\"Particles\">" << std::endl;
  fout.close();
}

template<typename T, typename FUNCTOR, vtkType VTKTYPE>
void VTKwriter<T,FUNCTOR,VTKTYPE>::closeVTU(
  const std::string& fullNamePiece)
{
  std::ofstream fout(fullNamePiece.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePiece << std::endl;
  }
  fout << "</Piece>" << std::endl;
  fout << "</UnstructuredGrid>\n";
  fout << "</VTKFile>\n";
  fout.close();
}




//TODO: make dimension insensitive
template<typename T, typename FUNCTOR, vtkType VTKTYPE>
void VTKwriter<T,FUNCTOR,VTKTYPE>::preambleVTI (
    const std::string& fullName,
    const LatticeR<D> extent0, const LatticeR<D> extent1,
    PhysR<T,D> origin, T delta)
{
  std::ofstream fout(fullName, std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }
  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"ImageData\" version=\"0.1\" ";
  if (_compress) {
    fout << "byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
  }
  else {
    fout << "byte_order=\"LittleEndian\">\n";
  }

  fout << "<ImageData WholeExtent=\""
       << extent0[0] <<" "<< extent1[0];
  for (int iDim=1; iDim<D; ++iDim){
    fout << " " << extent0[iDim] << " " << extent1[iDim];
  }
  fout << "\" Origin=\"" << origin[0];
  for (int iDim=1; iDim<D; ++iDim){
    fout << " " << origin[iDim];
  }
  fout << "\" Spacing=\"" << delta << " " << delta << " " << delta << "\">\n";
  fout << "<Piece Extent=\""
       << extent0[0] <<" "<< extent1[0];
  for (int iDim=1; iDim<D; ++iDim){
    fout << " " << extent0[iDim] <<" "<< extent1[iDim];
  }
  fout <<"\">\n";
  fout << "<PointData>\n";
  fout.close();
}


template<typename T, typename FUNCTOR, vtkType VTKTYPE>
void VTKwriter<T,FUNCTOR,VTKTYPE>::closeVTI(
  const std::string& fullNamePiece)
{
  std::ofstream fout(fullNamePiece, std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fullNamePiece << std::endl;
  }
  fout << "</ImageData>\n";
  fout << "</VTKFile>\n";
  fout.close();
}

template<typename T, typename FUNCTOR, vtkType VTKTYPE>
void VTKwriter<T,FUNCTOR,VTKTYPE>::closePiece(const std::string& fullNamePiece)
{
  std::ofstream fout(fullNamePiece, std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fullNamePiece << std::endl;
  }
  fout << "</PointData>\n";
  fout << "</Piece>\n";
  fout.close();
}


template<typename T, typename FUNCTOR, vtkType VTKTYPE>
void VTKwriter<T,FUNCTOR,VTKTYPE>::preamblePVD(
  const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }

  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"Collection\" version=\"0.1\" "
       << "byte_order=\"LittleEndian\">\n" << "<Collection>\n";
  fout.close();
}

template<typename T, typename FUNCTOR, vtkType VTKTYPE>
void VTKwriter<T,FUNCTOR,VTKTYPE>::closePVD(
  const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "</Collection>\n";
  fout << "</VTKFile>\n";
  fout.close();
}




template<typename T, typename FUNCTOR, vtkType VTKTYPE>
void VTKwriter<T,FUNCTOR,VTKTYPE>::preambleVTM(const std::string& fullNameVTM)
{
  std::ofstream fout(fullNameVTM, std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullNameVTM << std::endl;
  }
  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" "
       << "byte_order=\"LittleEndian\">\n"
       << "<vtkMultiBlockDataSet>\n" ;
  fout.close();
}

template<typename T, typename FUNCTOR, vtkType VTKTYPE>
void VTKwriter<T,FUNCTOR,VTKTYPE>::closeVTM(const std::string& fullNameVTM)
{
  std::ofstream fout(fullNameVTM, std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fullNameVTM << std::endl;
  }
  fout << "</vtkMultiBlockDataSet>\n";
  fout << "</VTKFile>\n";
  fout.close();
}


template<typename T, typename FUNCTOR, vtkType VTKTYPE>
void VTKwriter<T,FUNCTOR,VTKTYPE>::dataVTM(int iC, const std::string& fullNameVTM,
                                    const std::string& namePiece)
{
  std::ofstream fout(fullNameVTM, std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNameVTM << std::endl;
  }
  fout << "<Block index=\"" << iC << "\" >\n";
  fout << "<DataSet index= \"0\" " << "file=\"" << namePiece << "\">\n"
       << "</DataSet>\n";
  fout << "</Block>\n";
  fout.close();
}

template<typename T, typename FUNCTOR, vtkType VTKTYPE>
void VTKwriter<T,FUNCTOR,VTKTYPE>::dataPVD(int iT,
    const std::string& fullNamePVD, const std::string& namePiece)
{
  std::ofstream fout(fullNamePVD.c_str(),
                     std::ios::in | std::ios::out | std::ios::ate);
  if (fout) {
    fout.seekp(-25, std::ios::end); // jump -25 form the end of file to overwrite closePVD

    fout << "<DataSet timestep=\"" << iT << "\" "
         << "group=\"\" part=\"\" "
        << "file=\"" << namePiece << "\"/>\n";
    fout.close();
    closePVD(fullNamePVD);
  }
  else {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
}

template<typename T, typename FUNCTOR, vtkType VTKTYPE>
template<unsigned sourceDim>
void VTKwriter<T,FUNCTOR,VTKTYPE>::dataArrayPoints(
  const std::string& fullName,
  Vector<int,sourceDim> extent1, int iC )

{
  std::ofstream fout(fullName.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }
  fout << "<Points>" << std::endl;
  // Call dataArray with first functor (which should contain points)
  this->dataArraySingleFunctor(fullName, *_pointerVec[0], extent1, iC);
  fout << "</Points>" << std::endl;
}



template<typename T, typename FUNCTOR, vtkType VTKTYPE>
void VTKwriter<T,FUNCTOR,VTKTYPE>::cellDataVTU( const std::string& fullName, Vector<int, 1> extent1 ){
  std::ofstream fout(fullName.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }
  fout << "</PointData>" << std::endl;
  fout << "<CellData /> " << std::endl; //TODO: open tag missing?
  fout << "<Cells>" << std::endl;
  fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
       << std::endl;
  int32_t i32 = 0;
  for (int iTmp=0; iTmp<extent1[0]; ++iTmp){ fout << i32++ << " "; }
  fout << "</DataArray>" << std::endl;
  fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
       << std::endl;
  i32 = 1;
  for (int iTmp=0; iTmp<extent1[0]; ++iTmp){ fout << i32++ << " "; }
  fout << "</DataArray>" << std::endl;
  fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
       << std::endl;
  for (int iTmp=0; iTmp<extent1[0]; ++iTmp){ fout << 1 << " "; }
  fout << "</DataArray>" << std::endl;
  fout << "</Cells>" << std::endl;
}



template<typename T, typename FUNCTOR, vtkType VTKTYPE>
template<unsigned sourceDim>
void VTKwriter<T,FUNCTOR,VTKTYPE>::dataArraySingleFunctor(
  const std::string& fullName,
  FUNCTOR& f,
  Vector<int,sourceDim> extent1, int iC)
{
  std::ofstream fout(fullName.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }
  fout << "<DataArray type=\"Float32\" Name=\"" << f.getName() << "\" NumberOfComponents=\"" << f.getTargetDim() << "\" ";
  if (_compress || _binary) {
    fout << "format=\"binary\" encoding=\"base64\">\n";
  }
  else {
    fout << ">\n";
  }

  // Create functor input and output
  int i[4] = {0}; //4 as maximum dimension
  std::vector<T> evaluated(f.getTargetDim());

  //Set ic if parallel
  if constexpr(parallel){ i[0] = iC;}
  
  //Create float array
  size_t numberOfFloats;
  std::unique_ptr<float[]> streamFloat;
  int itter = 0;

  //VTI
  if constexpr(VTKTYPE==VTI){
    //2D
    if constexpr(D==2){
      numberOfFloats = f.getTargetDim() * (extent1[0]+2) * (extent1[1]+2);
      streamFloat = std::make_unique<float[]>(numberOfFloats);
      for (i[2]=-1; i[2]<extent1[1]+1; ++i[2]) {
      for (i[1]=-1; i[1]<extent1[0]+1; ++i[1]) {
        f(evaluated.data(),i);
        for (int iDim=0; iDim < f.getTargetDim(); ++iDim){
          streamFloat[itter] = float( evaluated[iDim] );
          ++itter;
        }
      }
      }
    //3D
    } else if constexpr(D==3){
      numberOfFloats = f.getTargetDim() * (extent1[0]+2) * (extent1[1]+2) * (extent1[2]+2);
      streamFloat = std::make_unique<float[]>(numberOfFloats);
      for (i[3]=-1; i[3]<extent1[2]+1; ++i[3]) {
      for (i[2]=-1; i[2]<extent1[1]+1; ++i[2]) {
      for (i[1]=-1; i[1]<extent1[0]+1; ++i[1]) {
        f(evaluated.data(),i);
        for (int iDim=0; iDim < f.getTargetDim(); ++iDim){
          streamFloat[itter] = float( evaluated[iDim] );
          ++itter;
        }
      }
      }
      }
    } else {
      clout << "Error: only 2D and 3D supportet" << std::endl;
    }

  //VTU
  } else if constexpr(VTKTYPE==VTU){
    numberOfFloats = f.getTargetDim() * extent1[0];
    streamFloat = std::make_unique<float[]>(numberOfFloats);
    //PARALLEL
    if constexpr(parallel){
      for (i[1]=0; i[1]<extent1[0]; ++i[1]) {
        f(evaluated.data(),i);
        for (int iDim=0; iDim < f.getTargetDim(); ++iDim){
          streamFloat[itter] = float( evaluated[iDim] );
          ++itter;
        }
      }
    //NONPARALLEL
    } else {
      for (i[0]=0; i[0]<extent1[0]; ++i[0]) {
        f(evaluated.data(),i);
        for (int iDim=0; iDim < f.getTargetDim(); ++iDim){
          streamFloat[itter] = float( evaluated[iDim] );
          ++itter;
        }
      }
    }
  } else {
    clout << "Error: VTP format not implemented yet" << std::endl;
  }

  //Define binarSize from number of floats
  uint32_t binarySize = static_cast<uint32_t>( numberOfFloats*sizeof(float) );

  //Write data (TODO: repair for non parallised vtu version)
  if (_compress) {
    // char buffer for functor data
    const unsigned char* charData = reinterpret_cast<unsigned char*>(streamFloat.get());
    // buffer for compression
    std::unique_ptr<unsigned char[]> comprData(new unsigned char[ binarySize ]);    // stack may be too small
 
    // compress data (not yet decoded as base64) by zlib
    uLongf sizeCompr = compressBound(binarySize);
    compress2( comprData.get(), &sizeCompr, charData, binarySize, -1);
 
    // encode prefix to base64 documented in  http://www.earthmodels.org/software/vtk-and-paraview/vtk-file-formats
    Base64Encoder<uint32_t> prefixEncoder(fout, 4);
    uint32_t prefix[4] = {1,binarySize,binarySize,static_cast<uint32_t>(sizeCompr)};
    prefixEncoder.encode(prefix, 4);
 
    // encode compressed data to base64
    Base64Encoder<unsigned char> dataEncoder( fout, sizeCompr );
    dataEncoder.encode(comprData.get(), sizeCompr);
  } else if (_binary) {
    // encode prefix to base64 documented in  http://www.earthmodels.org/software/vtk-and-paraview/vtk-file-formats
    Base64Encoder<uint32_t> prefixEncoder(fout, 1);
    prefixEncoder.encode(&binarySize, 1);
    //  write numbers from functor
    Base64Encoder<float> dataEncoder(fout, numberOfFloats);
    dataEncoder.encode(streamFloat.get(),numberOfFloats);
  } else {
    for ( size_t iOut = 0; iOut < numberOfFloats; ++iOut ) {
      fout << streamFloat[iOut] << " ";
    }
  }
  fout.close();
 
  std::ofstream ffout( fullName,  std::ios::out | std::ios::app );
  if (!ffout) {
    clout << "Error: could not open " << fullName << std::endl;
  }
  ffout << "\n</DataArray>\n";
  ffout.close();
}


}  // namespace OLB

#endif
