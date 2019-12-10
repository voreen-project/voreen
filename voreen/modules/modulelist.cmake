####################################################################################################
#  Module activation options                                                                       #   
####################################################################################################

# standard modules without external dependencies (enabled by default)
OPTION(VRN_MODULE_BASE                              "Provides basic rendering and data processing functionality"                                ON )
OPTION(VRN_MODULE_CONNEXE                           "Connected component analysis using the Connexe library"                                    ON )
OPTION(VRN_MODULE_EXPERIMENTAL                      "Provides experimental code not yet ready for productive use"                               OFF)
OPTION(VRN_MODULE_FLOWREEN                          "Provides flow visualization techniques"                                                    ON )
OPTION(VRN_MODULE_PLOTTING                          "Provides plotting components"                                                              ON )
OPTION(VRN_MODULE_POI                               "Points of interest in 2D and 3D"                                                           OFF)
OPTION(VRN_MODULE_PVM                               "Provides a volume reader for Stefan Roettger's PVM format"                                 ON )
OPTION(VRN_MODULE_SAMPLE                            "Contains example processors"                                                               OFF)
OPTION(VRN_MODULE_SEGY                              "Provides a volume reader for the SEGY format"                                              ON )
OPTION(VRN_MODULE_STEREOSCOPY                       "Provides stereo rendering"                                                                 OFF)
OPTION(VRN_MODULE_STAGING                           "Experimental code almost ready for productive use"                                         ON )
OPTION(VRN_MODULE_SURFACE                           "Provides isosurface extraction"                                                            ON )
OPTION(VRN_MODULE_RANDOMWALKER                      "3D random walker implementation"                                                           ON )
OPTION(VRN_MODULE_VOLUMELABELING                    "Provides interactive volume labeling (requires OpenGL compatiblity profile)"               OFF )

# important modules with external dependencies that are usually available on a target system
OPTION(VRN_MODULE_DEVIL                             "Provides image I/O and processing using the DevIL library"                                 ON )
OPTION(VRN_MODULE_ZIP                               "Provides ZIP support via the zlib library"                                                 ON )
OPTION(VRN_MODULE_PYTHON                            "Python bindings, allowing scripts to access the Voreen API"                                ON )

# modules with non-standard external dependencies
OPTION(VRN_MODULE_FFMPEG                            "Provides video encoding using the FFmpeg library"                                          OFF)
OPTION(VRN_MODULE_HDF5                              "Provides support for reading and writing HDF5 format volumes"                              OFF)
OPTION(VRN_MODULE_TIFF                              "Volume readers for multi-image TIFF files and OME-TIFF stacks"                             OFF)
OPTION(VRN_MODULE_OPENCL                            "Provides an OpenCL wrapper and OpenCL-based processors"                                    OFF)
OPTION(VRN_MODULE_OPENMP                            "Activates OpenMP code for parallelization"                                                 OFF)
OPTION(VRN_MODULE_BIGDATAIMAGEPROCESSING            "Provides image processing support for very large volumes"                                  OFF)
OPTION(VRN_MODULE_VESSELNETWORKANALYSIS             "Provides processors for analyzing volumetric vessel structures"                            OFF)

# dicom support
OPTION(VRN_MODULE_GDCM                              "Provides a DICOM volume reader using the GDCMv2 library"                                   OFF)

# itk support
OPTION(VRN_MODULE_ITK                               "Wrapper for ITK filters as Voreen processors"                                              OFF)
OPTION(VRN_MODULE_ITK_GENERATED                     "Automatically generated Voreen processors as wrappers for ITK filters"                     OFF)

# additional modules
OPTION(VRN_MODULE_DEPRECATED                        "Deprecated code to be removed in future releases"                                          OFF)
