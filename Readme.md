# Voreen - The Volume Rendering Engine

Voreen is an open source rapid application development framework for the interactive visualization 
and analysis of multi-modal volumetric data sets. It provides GPU-based volume rendering and data 
analysis techniques and offers high flexibility when developing new analysis workflows in 
collaboration with domain experts. The Voreen framework consists of a multi-platform C++ library, 
which can be easily integrated into existing applications, and a Qt-based stand-alone application. 
It is licensed under the terms of the GNU General Public License version 2.

The Voreen project has been initiated and was originally maintained by the Visualization & Computer 
Graphics Research Group at the University of Münster as part of the collaborative research center 
SFB 656 MoBil (Project Z1, Project Ö). Since 2018, Voreen is collaboratively developed by the Pattern 
Recognition and Image Analysis (PRIA) Research Group and the VISualization & graphIX (VISIX) Research Group.

For build instructions and further information, please refer to the [project website](http://voreen.uni-muenster.de).

## Directory Structure


| directory | purpose |
| --------- | ------- |
|include                         | framework headers |
|src                             | framework source files |
|modules                         | standard plugin modules |
|custommodules                   | location for custom (third-party) modules |
|bin                             | location where compiled binaries are written to |
| | |
|apps                            | applications using voreen |
|apps/simple                     | simple Qt- and GLUT-based applications |
|apps/voreenve                   | main Voreen application - the Voreen Visualization Environment |
|apps/voreentool                 | command-line interface to the Voreen library, executing workspaces/networks |
| | |
|data                            | location for user-generated data (log files, settings, screenshots, ...) |
| | |
|resource                        | contains read-only resources, grouped by the corresponding Voreen component (i.e., voreencore, voreenqt, voreenve) |
|resource/<...>/fonts            | (possibly) necessary fonts |
|resource/<...>/geometry         | sample geometry data |
|resource/<...>/scripts          | Python scripts which can be executed from within Voreen |
|resource/<...>/transferfuncs    | pre-defined transfer functions |
|resource/<...>/volumes          | sample volume data |
|resource/<...>/workspaces       | sample Voreen workspaces |
| | |
|doc                             | documentation |
|ext                             | external dependencies to other libraries |
