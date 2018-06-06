VoreenVTK Build Instructions
----------------------------

* Install Cmake (www.cmake.org) and VTK (www.vtk.org)
  - Guide: http://vtkblog.blogspot.com/2007/04/install-build-vtk-from-source-in-visual.html
  - Make sure to build VTK as shared library (cmake option)

* Compile Voreen core library

* Run Cmake on this directory (apps/voreenvtk):
  - press 'configure'
  - specify Voreen root directory
  - select external libraries (must match the configuration the Voreen lib has been built with)
  - press 'generate'

* Open generated VS solution file (.sln) and compile
  - If 'error LNK2005' occurs (multiple definitions), VTK has most likely not been built as shared lib
    -> Rebuild VTK, or (NOT RECOMMENDED!): add /FORCE:MULTIPLE to linker options of 'TestProject' application

* Copy VTK dlls to build directory

* The workspace to load can be passed as command line argument
