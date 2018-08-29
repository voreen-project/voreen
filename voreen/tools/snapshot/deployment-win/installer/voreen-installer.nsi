;
; Voreen NSIS installer script (MUI version). 
;
; Installs the application into a user-specified directory, generates an uninstaller as well as entries in the start menu 
; and in "Add/Remove Programs". The installer does not require administrator privileges.
;
; Customize the script by setting up the "Configuration" section and adding the files to be installed to the "Files" section.
; No further modifications to this file should be necessary.
;
; For compilation place this script in the application's root directory and make sure that all Voreen files 
; as well as the installer bitmaps specified in the "Interface Settings" section are available.
;

;--------------------------------
; Configuration

  !define AppName "VoreenVE"               ;  displayed in the installer and used as app directory name and shortcut name       
  !define Version "5.0.1"                  ;  added to appname    
  !define ExecutableName "voreenve.exe"    ;  filename of the executable
  
  ; settings for the "Add/Remove Programs" entry
  !define DisplayName "VoreenVE ${Version}– Visualization Environment"
  !define Publisher "University of Muenster, Department of Computer Science"
  !define HelpLink "https://www.uni-muenster.de/Voreen/"

;--------------------------------
; Include Modern UI

  !include "MUI2.nsh"
  
;--------------------------------
; Basic NSIS settings
  
  ;Name and file
  Name "${AppName} ${Version}"
  OutFile "${AppName}-${Version}-Installer.exe"

  ;Default installation folder
  InstallDir "$PROGRAMFILES\${AppName}-${Version}"
  
  ;Request application privileges for Windows Vista
  RequestExecutionLevel user

;--------------------------------
; Interface Settings

  !define MUI_ABORTWARNING
  !define MUI_HEADERIMAGE
  !define MUI_HEADERIMAGE_BITMAP "installer-header.bmp"
  !define MUI_WELCOMEFINISHPAGE_BITMAP "wizard-blue.bmp"
  !define MUI_UNWELCOMEFINISHPAGE_BITMAP "wizard-uninstall-blue.bmp"

;--------------------------------
; Installer Pages (which pages are to be shown during the installation/uninstallation workflow)

  !insertmacro MUI_PAGE_WELCOME
  !insertmacro MUI_PAGE_LICENSE "..\..\..\..\LICENSE-academic.txt"
  #!insertmacro MUI_PAGE_COMPONENTS
  !insertmacro MUI_PAGE_DIRECTORY
  !insertmacro MUI_PAGE_INSTFILES
  !define MUI_FINISHPAGE_SHOWREADME
  !define MUI_FINISHPAGE_SHOWREADME_TEXT "Create Shortcut on Desktop"
  !define MUI_FINISHPAGE_SHOWREADME_FUNCTION CreateDesktopShortcut
  !define MUI_FINISHPAGE_RUN "$INSTDIR\${ExecutableName}"
  !insertmacro MUI_PAGE_FINISH
  
  !insertmacro MUI_UNPAGE_WELCOME
  !insertmacro MUI_UNPAGE_CONFIRM
  !insertmacro MUI_UNPAGE_INSTFILES
  !insertmacro MUI_UNPAGE_FINISH

;--------------------------------
; Languages
 
  !insertmacro MUI_LANGUAGE "English"

;--------------------------------
; Macros for automated uninstallation of registered files (DO NOT TOUCH!)
  
!define UninstLog "uninstall.log"
Var UninstLog

; Uninstall log file missing.
LangString UninstLogMissing ${LANG_ENGLISH} "${UninstLog} not found!$\r$\nUninstallation cannot proceed!"
 
; AddItem macro
!macro AddItem Path
 FileWrite $UninstLog "${Path}$\r$\n"
!macroend
!define AddItem "!insertmacro AddItem"
 
; File macro
!macro File FilePath FileName
 IfFileExists "$OUTDIR\${FileName}" +2
  FileWrite $UninstLog "$OUTDIR\${FileName}$\r$\n"
 File "${FilePath}${FileName}"
!macroend
!define File "!insertmacro File"
 
; Copy files macro
!macro CopyFiles SourcePath DestPath
 IfFileExists "${DestPath}" +2
  FileWrite $UninstLog "${DestPath}$\r$\n"
 CopyFiles "${SourcePath}" "${DestPath}"
!macroend
!define CopyFiles "!insertmacro CopyFiles"
 
; Rename macro
!macro Rename SourcePath DestPath
 IfFileExists "${DestPath}" +2
  FileWrite $UninstLog "${DestPath}$\r$\n"
 Rename "${SourcePath}" "${DestPath}"
!macroend
!define Rename "!insertmacro Rename"
 
; CreateDirectory macro
!macro CreateDirectory Path
 CreateDirectory "${Path}"
 FileWrite $UninstLog "${Path}$\r$\n"
!macroend
!define CreateDirectory "!insertmacro CreateDirectory"
 
; SetOutPath macro
!macro SetOutPath Path
 SetOutPath "${Path}"
 FileWrite $UninstLog "${Path}$\r$\n"
!macroend
!define SetOutPath "!insertmacro SetOutPath"
 
; WriteUninstaller macro
!macro WriteUninstaller Path
 WriteUninstaller "${Path}"
 FileWrite $UninstLog "${Path}$\r$\n"
!macroend
!define WriteUninstaller "!insertmacro WriteUninstaller"
 
Section -openlogfile
 CreateDirectory "$INSTDIR"
 IfFileExists "$INSTDIR\${UninstLog}" +3
  FileOpen $UninstLog "$INSTDIR\${UninstLog}" w
 Goto +4
  SetFileAttributes "$INSTDIR\${UninstLog}" NORMAL
  FileOpen $UninstLog "$INSTDIR\${UninstLog}" a
  FileSeek $UninstLog 0 END
SectionEnd

; Attempt to give the UAC plug-in a user process and an admin process.
Function .OnInit
 
UAC_Elevate:
    UAC::RunElevated 
    StrCmp 1223 $0 UAC_ElevationAborted ; UAC dialog aborted by user?
    StrCmp 0 $0 0 UAC_Err ; Error?
    StrCmp 1 $1 0 UAC_Success ;Are we the real deal or just the wrapper?
    Quit
 
UAC_Err:
    MessageBox mb_iconstop "Unable to elevate, error $0"
    Abort
 
UAC_ElevationAborted:
    # elevation was aborted, run as normal?
    MessageBox mb_iconstop "This installer requires admin access, aborting!"
    Abort
 
UAC_Success:
    StrCmp 1 $3 +4 ;Admin?
    StrCmp 3 $1 0 UAC_ElevationAborted ;Try again?
    MessageBox mb_iconstop "This installer requires admin access, try again"
    goto UAC_Elevate 
 
FunctionEnd

Function .OnInstFailed
    UAC::Unload ;Must call unload!
FunctionEnd
 
Function .OnInstSuccess
    UAC::Unload ;Must call unload!
FunctionEnd


; Attempt to give the UAC plug-in a user process and an admin process. (Uninstaller)
Function un.OnInit
 
UAC_Elevate:
    UAC::RunElevated 
    StrCmp 1223 $0 UAC_ElevationAborted ; UAC dialog aborted by user?
    StrCmp 0 $0 0 UAC_Err ; Error?
    StrCmp 1 $1 0 UAC_Success ;Are we the real deal or just the wrapper?
    Quit
 
UAC_Err:
    MessageBox mb_iconstop "Unable to elevate, error $0"
    Abort
 
UAC_ElevationAborted:
    # elevation was aborted, run as normal?
    MessageBox mb_iconstop "This installer requires admin access, aborting!"
    Abort
 
UAC_Success:
    StrCmp 1 $3 +4 ;Admin?
    StrCmp 3 $1 0 UAC_ElevationAborted ;Try again?
    MessageBox mb_iconstop "This installer requires admin access, try again"
    goto UAC_Elevate 
 
FunctionEnd

; --------------------------------------------------------------------------------
; Files

Section "Main Section" SecMain

    !include voreen-installer-files.nsi
    
  ; ${SetOutPath} "$INSTDIR"

  ; ### ADD INSTALLATION FILES HERE ###
  
  ; ; App files (root dir)
  ; ${File} "" ${ExecutableName}
; ;  ${File} "" "voltool.exe"
  ; ${File} "" "LICENSE.txt"
  ; ${File} "" "LICENSE-academic.txt"
  ; ${File} "" "Changelog.txt"
  ; ${File} "" "CREDITS.txt"
  
  ; ; Qt dlls
  ; ${File} "" "QtCore4.dll"
  ; ${File} "" "QtGui4.dll"
  ; ${File} "" "QtOpenGL4.dll"
  
  ; ; VS 2008 Redistributables
  ; ${File} "" "msvcm90.dll"
  ; ${File} "" "msvcp90.dll"
  ; ${File} "" "msvcr90.dll"
  ; ${File} "" "Microsoft.VC90.CRT.manifest"
  
  ; ; Libraries
  ; ${File} "" "avcodec-52.dll"
  ; ${File} "" "avdevice-52.dll"
  ; ${File} "" "avformat-52.dll"
  ; ${File} "" "avutil-50.dll"
  ; ${File} "" "swscale-0.dll"
  ; ${File} "" "python26.dll"
  ; ${File} "" "DevIL.dll"
  ; ${File} "" "ILU.dll"
  ; ${File} "" "freetype6.dll"
  ; ${File} "" "ftgl.dll"
  ; ${File} "" "jpeg62.dll"
  ; ${File} "" "libtiff3.dll"
  ; ${File} "" "zlib1.dll"
    
  ; ; directory data/
  ; ${AddItem} "$INSTDIR\data"
  ; ${SetOutPath} "$INSTDIR\data"

  ; ; directory data/fonts
  ; ${AddItem} "$INSTDIR\data\fonts"
  ; ${SetOutPath} "$INSTDIR\data\fonts"
  ; ${File} "data\fonts\" "Vera.ttf"
  ; ${File} "data\fonts\" "VeraMono.ttf"
  
  ; ; directory data/plotting
  ; ${AddItem} "$INSTDIR\data\plotting"
  ; ${SetOutPath} "$INSTDIR\data\plotting"
  ; ${File} "data\plotting\" "data_lineplot.csv"
  
  ; ; directory data/networks
  ; ${AddItem} "$INSTDIR\data\networks"
  ; ${SetOutPath} "$INSTDIR\data\networks"
    
  ; ; directory data/scripts
  ; ${AddItem} "$INSTDIR\data\scripts"
  ; ${SetOutPath} "$INSTDIR\data\scripts"
  ; ${File} "data\scripts\" "fps.py"
  ; ${File} "data\scripts\" "fps-clock.py"
  ; ${File} "data\scripts\" "snapshot.py"
  
  ; ; directory data/textures
  ; ${AddItem} "$INSTDIR\data\textures"
  ; ${SetOutPath} "$INSTDIR\data\textures"
  ; ${File} "data\textures\" "error.tga"
  ; ${File} "data\textures\" "anterior2.png"
  ; ${File} "data\textures\" "axial_b.png"
  ; ${File} "data\textures\" "axial_t.png"
  ; ${File} "data\textures\" "coronal_b.png"
  ; ${File} "data\textures\" "coronal_f.png"
  ; ${File} "data\textures\" "inferior2.png"
  ; ${File} "data\textures\" "lateral2.png"
  ; ${File} "data\textures\" "sagittal_l.png"
  ; ${File} "data\textures\" "sagittal_r.png"
  ; ${File} "data\textures\" "septal2.png"
  
  ; ; directory data/transferfuncs
  ; ${AddItem} "$INSTDIR\data\transferfuncs"
  ; ${SetOutPath} "$INSTDIR\data\transferfuncs"
  ; ${File} "data\transferfuncs\" "nucleon.tfi"
  ; ${File} "data\transferfuncs\" "nucleon2.tfi"
  ; ${File} "data\transferfuncs\" "walnut.tfi"
  ; ${File} "data\transferfuncs\" "chromadepthspectrum.bmp"
  
  ; ; directory data/volumes
  ; ${AddItem} "$INSTDIR\data\volumes"
  ; ${SetOutPath} "$INSTDIR\data\volumes"
  ; ${File} "data\volumes\" "nucleon.dat"
  ; ${File} "data\volumes\" "nucleon.raw"
  ; ${File} "data\volumes\" "walnut.dat"
  ; ${File} "data\volumes\" "walnut.raw"
  ; ${File} "data\volumes\" "walnut-transformed.dat"
  ; ${File} "data\volumes\" "walnut-transformed.raw"
  
  ; ; directory data/workspaces
  ; ${AddItem} "$INSTDIR\data\workspaces"
  ; ${SetOutPath} "$INSTDIR\data\workspaces"
  ; ${File} "data\workspaces\" "animation-simple.vws"
  ; ${File} "data\workspaces\" "arbitraryclipping.vws"
  ; ${File} "data\workspaces\" "explodedviews.vws"
  ; ${File} "data\workspaces\" "fancy.vws"
  ; ${File} "data\workspaces\" "glslraycaster.vws"
  ; ${File} "data\workspaces\" "multicanvas.vws"
  ; ${File} "data\workspaces\" "multivolume.vws"
  ; ${File} "data\workspaces\" "plotsample.vws"
  ; ${File} "data\workspaces\" "quadview.vws"
  ; ${File} "data\workspaces\" "raycasting-with-geometry.vws"
  ; ${File} "data\workspaces\" "renderloop.vws"
  ; ${File} "data\workspaces\" "simple.vws"
  ; ${File} "data\workspaces\" "sliceview.vws"
  ; ${File} "data\workspaces\" "standard.vws"
  ; ${File} "data\workspaces\" "walnut.vws"
  
  ; ; directory data/workspaces/tools
  ; ${AddItem} "$INSTDIR\data\workspaces\tools"
  ; ${SetOutPath} "$INSTDIR\data\workspaces\tools"
  ; ${File} "data\workspaces\tools\" "volumecropping.vws"  
  ; ${File} "data\workspaces\tools\" "volumeregistration.vws"  
  
  ; ; directory doc/
  ; ${AddItem} "$INSTDIR\doc"
  ; ${SetOutPath} "$INSTDIR\doc"
  ; ; directory doc/gettingstarted
  ; ${AddItem} "$INSTDIR\doc\gettingstarted"
  ; ${SetOutPath} "$INSTDIR\doc\gettingstarted"
  ; ${File} "doc\gettingstarted\" "gsg.html"
  ; ; directory doc/gettingstarted/images
  ; ${AddItem} "$INSTDIR\doc\gettingstarted\images"
  ; ${SetOutPath} "$INSTDIR\doc\gettingstarted\images"
  ; ${File} "doc\gettingstarted\images\" "image_numbers.png"
  ; ${File} "doc\gettingstarted\images\" "preview.png"
  ; ${File} "doc\gettingstarted\images\" "processor.png"
  ; ${File} "doc\gettingstarted\images\" "processors.png"
  ; ${File} "doc\gettingstarted\images\" "property.png"
  ; ${File} "doc\gettingstarted\images\" "tc.png"
  ; ${File} "doc\gettingstarted\images\" "toolbar.png"
  ; ${File} "doc\gettingstarted\images\" "tooltip.png"
  ; ${File} "doc\gettingstarted\images\" "vis_mode.png"
  ; ${File} "doc\gettingstarted\images\" "vss_proc.png"
  
   ; ; directory doc/animation
  ; ${AddItem} "$INSTDIR\doc\animation"
  ; ${SetOutPath} "$INSTDIR\doc\animation"
  ; ${File} "doc\animation\" "animation.html"
  ; ; directory doc/animation/images
  ; ${AddItem} "$INSTDIR\doc\animation\images"
  ; ${SetOutPath} "$INSTDIR\doc\animation\images"
  ; ${File} "doc\animation\images\" "overview.png"
  ; ${File} "doc\animation\images\" "voreen_cut2_200.png"
  ; ${File} "doc\animation\images\" "voreen_cut3_200.png"
  ; ${File} "doc\animation\images\" "voreen_cut_200.png"
  
  ; ; directory glsl/
  ; ; ${AddItem} "$INSTDIR\glsl"
  ; ; ${SetOutPath} "$INSTDIR\glsl"
  ; ; ${File} "glsl\" "colorcoding.vert"
  ; ; ${File} "glsl\" "colorcoding2d.frag"
  ; ; ${File} "glsl\" "colorcoding3d.frag" 
  ; ; ${File} "glsl\" "eep_clipping.frag"
  ; ; ${File} "glsl\" "eep_depth.frag"
  ; ; ${File} "glsl\" "eep_depth.vert"
  ; ; ${File} "glsl\" "eep_geometry.frag"
  ; ; ${File} "glsl\" "eep_inside_volume.frag"
  ; ; ${File} "glsl\" "eep_inside_volume.vert"
  ; ; ${File} "glsl\" "eep_jitter.frag"
  ; ; ${File} "glsl\" "eep_simple.frag"
  ; ; ${File} "glsl\" "eep_simple.vert"
  ; ; ${File} "glsl\" "eep_texcoord.vert"
  ; ; ${File} "glsl\" "eep_vertex.vert"  
  ; ; ${File} "glsl\" "gabor.frag"
  ; ; ${File} "glsl\" "mc_extract.vert"
  ; ; ${File} "glsl\" "mc_render.vert"
  ; ; ${File} "glsl\" "mc_render.frag"
  ; ; ${File} "glsl\" "mod_colorcoding.frag"
  ; ; ${File} "glsl\" "mod_phong.frag"
  ; ; ${File} "glsl\" "mod_phong.frag"
  ; ; ${File} "glsl\" "phong.frag"
  ; ; ${File} "glsl\" "phong.vert"
  ; ; ${File} "glsl\" "pp_background.frag"
  ; ; ${File} "glsl\" "pp_binary.frag"
  ; ; ${File} "glsl\" "pp_canny.frag"
  ; ; ${File} "glsl\" "pp_colordepth.frag"
  ; ; ${File} "glsl\" "pp_compositor.frag"
  ; ; ${File} "glsl\" "pp_convolution.frag"
  ; ; ${File} "glsl\" "pp_depthdarkening.frag"
  ; ; ${File} "glsl\" "pp_depthpeeling.frag"  
  ; ; ${File} "glsl\" "pp_depthpeeling.vert"  
  ; ; ${File} "glsl\" "pp_dilation.frag"  
  ; ; ${File} "glsl\" "pp_distance.frag"  
  ; ; ${File} "glsl\" "pp_edgedetect.frag"  
  ; ; ${File} "glsl\" "pp_erosion.frag"  
  ; ; ${File} "glsl\" "pp_fade.frag"  
  ; ; ${File} "glsl\" "pp_gaussian.frag"  
  ; ; ${File} "glsl\" "pp_grayscale.frag"  
  ; ; ${File} "glsl\" "pp_imagethreshold.frag"  
  ; ; ${File} "glsl\" "pp_labeling.frag"  
  ; ; ${File} "glsl\" "pp_mask.frag"  
  ; ; ${File} "glsl\" "pp_mean.frag"  
  ; ; ${File} "glsl\" "pp_median.frag"  
  ; ; ${File} "glsl\" "pp_nonminmax.frag"  
  ; ; ${File} "glsl\" "pp_orientationoverlay.frag"  
  ; ; ${File} "glsl\" "pp_pipethrough.frag"  
  ; ; ${File} "glsl\" "pp_scale.frag"  
  ; ; ${File} "glsl\" "pp_unary.frag"  
  ; ; ${File} "glsl\" "rc_curvature.frag"  
  ; ; ${File} "glsl\" "rc_firsthit.frag"  
  ; ; ${File} "glsl\" "rc_hitpoints.frag"  
  ; ; ${File} "glsl\" "rc_id.frag"  
  ; ; ${File} "glsl\" "rc_multivolume.frag"  
  ; ; ${File} "glsl\" "rc_rgb.frag"  
  ; ; ${File} "glsl\" "rc_segmentation.frag"    
  ; ; ${File} "glsl\" "rc_simple.frag"  
  ; ; ${File} "glsl\" "rc_singlevolume.frag"  
  ; ; ${File} "glsl\" "sl_base.frag"  
  ; ; ${File} "glsl\" "sl_base.vert"    
  ; ; ${File} "glsl\" "sl_frontcompositing.frag"  
  ; ; ${File} "glsl\" "sl_halfslicing.frag"  
  ; ; ${File} "glsl\" "sl_singlevolume.frag"    
  ; ; ${File} "glsl\" "sl_singlevolume.vert"  
  ; ; ${File} "glsl\" "spotnoise.frag"  
  ; ; ${File} "glsl\" "spotnoise2d.vert"  
  ; ; ${File} "glsl\" "spotnoise3d.vert"  
  ; ; ${File} "glsl\" "streamlinerenderer3d.frag"  
  ; ; ${File} "glsl\" "textoverlay.frag"    
  
  ; ; directory glsl/modules
  ; ; ${AddItem} "$INSTDIR\glsl\modules"
  ; ; ${SetOutPath} "$INSTDIR\glsl\modules"
  ; ; ${File} "glsl\modules\" "mod_bricking.frag"    
  ; ; ${File} "glsl\modules\" "mod_compositing.frag"    
  ; ; ${File} "glsl\modules\" "mod_curvature.frag"    
  ; ; ${File} "glsl\modules\" "mod_depth.frag"    
  ; ; ${File} "glsl\modules\" "mod_filtering.frag"    
  ; ; ${File} "glsl\modules\" "mod_firsthit.frag"    
  ; ; ${File} "glsl\modules\" "mod_gradients.frag"    
  ; ; ${File} "glsl\modules\" "mod_lightprobe.frag"    
  ; ; ${File} "glsl\modules\" "mod_masking.frag"    
  ; ; ${File} "glsl\modules\" "mod_normdepth.frag"    
  ; ; ${File} "glsl\modules\" "mod_raysetup.frag"    
  ; ; ${File} "glsl\modules\" "mod_sampler2d.frag"    
  ; ; ${File} "glsl\modules\" "mod_sampler3d.frag"    
  ; ; ${File} "glsl\modules\" "mod_segmentation.frag"    
  ; ; ${File} "glsl\modules\" "mod_shading.frag"    
  ; ; ${File} "glsl\modules\" "mod_shadows.frag"    
  ; ; ${File} "glsl\modules\" "mod_sketch.frag"    
  ; ; ${File} "glsl\modules\" "mod_transfunc.frag"    
  ; ; ${File} "glsl\modules\" "vrn_shaderincludes.frag"    
  
  ; ; directory glsl/modules/bricking
  ; ; ${AddItem} "$INSTDIR\glsl\modules\bricking"
  ; ; ${SetOutPath} "$INSTDIR\glsl\modules\bricking"
  ; ; ${File} "glsl\modules\bricking\" "mod_adaptive_sampling.frag"    
  ; ; ${File} "glsl\modules\bricking\" "mod_basics.frag"    
  ; ; ${File} "glsl\modules\bricking\" "mod_bricking.frag"    
  ; ; ${File} "glsl\modules\bricking\" "mod_global_variables.frag"    
  ; ; ${File} "glsl\modules\bricking\" "mod_interpolation.frag"    
  ; ; ${File} "glsl\modules\bricking\" "mod_lookups.frag"    
  ; ; ${File} "glsl\modules\bricking\" "mod_math.frag"    
  ; ; ${File} "glsl\modules\bricking\" "mod_uniforms.frag"    
   
  ; ; ; directory glsl/qt/rendertargetviewer
  ; ; ${AddItem} "$INSTDIR\glsl\qt"
  ; ; ${SetOutPath} "$INSTDIR\glsl\qt"
  ; ; ${AddItem} "$INSTDIR\glsl\qt\rendertargetviewer"
  ; ; ${SetOutPath} "$INSTDIR\glsl\qt\rendertargetviewer"
  ; ; ${File} "glsl\qt\rendertargetviewer\" "color.frag"    
  ; ; ${File} "glsl\qt\rendertargetviewer\" "inversecolor.frag"    

  ; ; ; directory glsl/utils
  ; ; ${AddItem} "$INSTDIR\glsl\utils"
  ; ; ${SetOutPath} "$INSTDIR\glsl\utils"
  ; ; ${File} "glsl\utils\" "blendwithimage.frag"   
  ; ; ${File} "glsl\utils\" "copyimage.frag"     
  ; ; ${File} "glsl\utils\" "passthrough.vert"   
  ; ; ${File} "glsl\utils\" "vrn_rendertexture.frag"   
  ; ; shaders end
  
  ; ; directory licenses/
  ; ${AddItem} "$INSTDIR\licenses"
  ; ${SetOutPath} "$INSTDIR\licenses"
 
  ; ${AddItem} "$INSTDIR\licenses\connexe"
  ; ${SetOutPath} "$INSTDIR\licenses\connexe"
  ; ${File} "licenses\connexe\" "COPYRIGHT"

  ; ${AddItem} "$INSTDIR\licenses\dcmtk"
  ; ${SetOutPath} "$INSTDIR\licenses\dcmtk"
  ; ${File} "licenses\dcmtk\" "COPYRIGHT.txt"
  
  ; ${AddItem} "$INSTDIR\licenses\devil"
  ; ${SetOutPath} "$INSTDIR\licenses\devil"
  ; ${File} "licenses\devil\" "lgpl.txt"
  
  ; ; ${AddItem} "$INSTDIR\licenses\fboclass"
  ; ; ${SetOutPath} "$INSTDIR\licenses\fboclass"
  ; ; ${File} "licenses\fboclass\" "license.txt"
  
  ; ${AddItem} "$INSTDIR\licenses\freetype"
  ; ${SetOutPath} "$INSTDIR\licenses\freetype"
  ; ${File} "licenses\freetype\" "FTL.TXT"
  
  ; ${AddItem} "$INSTDIR\licenses\ftgl"
  ; ${SetOutPath} "$INSTDIR\licenses\ftgl"
  ; ${File} "licenses\ftgl\" "license.txt"
  
  ; ${AddItem} "$INSTDIR\licenses\glew"
  ; ${SetOutPath} "$INSTDIR\licenses\glew"
  ; ${File} "licenses\glew\" "license.txt"
  
  ; ${AddItem} "$INSTDIR\licenses\hpmc"
  ; ${SetOutPath} "$INSTDIR\licenses\hpmc"
  ; ${File} "licenses\libtiff\" "license.txt"
  
  ; ${AddItem} "$INSTDIR\licenses\libjpeg"
  ; ${SetOutPath} "$INSTDIR\licenses\libjpeg"
  ; ${File} "licenses\libjpeg\" "license.txt"
  
  ; ${AddItem} "$INSTDIR\licenses\libpng"
  ; ${SetOutPath} "$INSTDIR\licenses\libpng"
  ; ${File} "licenses\libpng\" "libpng-LICENSE.txt"
  
  ; ${AddItem} "$INSTDIR\licenses\libtiff"
  ; ${SetOutPath} "$INSTDIR\licenses\libtiff"
  ; ${File} "licenses\libtiff\" "license.txt"

  ; ${AddItem} "$INSTDIR\licenses\python"
  ; ${SetOutPath} "$INSTDIR\licenses\python"
  ; ${File} "licenses\python\" "LICENSE"

  ; ${AddItem} "$INSTDIR\licenses\Qt"
  ; ${SetOutPath} "$INSTDIR\licenses\Qt"
  ; ${File} "licenses\Qt\" "LICENSE"
  
  ; ${AddItem} "$INSTDIR\licenses\tinyxml"
  ; ${SetOutPath} "$INSTDIR\licenses\tinyxml"
  ; ${File} "licenses\tinyxml\" "license.txt"
  ;
  ; ${AddItem} "$INSTDIR\licenses\tinyobj"
  ; ${SetOutPath} "$INSTDIR\licenses\tinyobj"
  ; ${File} "licenses\tinyobj\" "license.txt"

  ; ${AddItem} "$INSTDIR\licenses\v3"
  ; ${SetOutPath} "$INSTDIR\licenses\v3"
  ; ${File} "licenses\v3\" "LICENSE.txt"
 
  ; ${AddItem} "$INSTDIR\licenses\zlib"
  ; ${SetOutPath} "$INSTDIR\licenses\zlib"
  ; ${File} "licenses\zlib\" "license.txt"
  ; ; licenses end
  
  ### INSTALLATION FILES END ###

  ; Create uninstaller
  WriteUninstaller "$INSTDIR\Uninstall.exe"

  ; Startmenu entries
  SetOutPath "$INSTDIR"
  CreateDirectory "$SMPROGRAMS\${AppName} ${Version}"
  createShortCut "$SMPROGRAMS\${AppName} ${Version}\${AppName} ${Version}.lnk" "$INSTDIR\${ExecutableName}"
  createShortCut "$SMPROGRAMS\${AppName} ${Version}\Uninstall.lnk" "$INSTDIR\Uninstall.exe"
  
  ; Add entry to "Add/Remove Programs"
  WriteRegStr HKCU  "Software\Microsoft\Windows\CurrentVersion\Uninstall\${AppName}-${Version}" \
                 "DisplayName" "${DisplayName}"
  WriteRegStr HKCU  "Software\Microsoft\Windows\CurrentVersion\Uninstall\${AppName}-${Version}" \
                 "UninstallString" "$\"$INSTDIR\Uninstall.exe$\""
  WriteRegStr HKCU  "Software\Microsoft\Windows\CurrentVersion\Uninstall\${AppName}-${Version}" \
                 "Publisher" "${Publisher}"
  WriteRegStr HKCU  "Software\Microsoft\Windows\CurrentVersion\Uninstall\${AppName}-${Version}" \
                 "HelpLink" "${HelpLink}"
 
SectionEnd
  
Function CreateDesktopShortcut
   createShortcut "$DESKTOP\${AppName} ${Version}.lnk" "$INSTDIR\${ExecutableName}"
FunctionEnd

;--------------------------------
; Uninstaller Section (DO NOT TOUCH!)
  
Section -closelogfile
 FileClose $UninstLog
 SetFileAttributes "$INSTDIR\${UninstLog}" READONLY|SYSTEM|HIDDEN
SectionEnd
 
Section Uninstall
 
 ; first remove start menu entries
 Delete  "$SMPROGRAMS\${AppName} ${Version}\${AppName} ${Version}.lnk"
 Delete  "$SMPROGRAMS\${AppName} ${Version}\Uninstall.lnk"
 RMDir "$SMPROGRAMS\${AppName} ${Version}"
 Delete  "$DESKTOP\${AppName} ${Version}.lnk"
 
 ; delete registry entry for "Add/Remove Programs" 
 DeleteRegKey HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${AppName} ${Version}"
 
 ; Can't uninstall if uninstall log is missing!
 IfFileExists "$INSTDIR\${UninstLog}" +3
  MessageBox MB_OK|MB_ICONSTOP "$(UninstLogMissing)"
   Abort
 
 Push $R0
 Push $R1
 Push $R2
 SetFileAttributes "$INSTDIR\${UninstLog}" NORMAL
 FileOpen $UninstLog "$INSTDIR\${UninstLog}" r
 StrCpy $R1 0
 
 GetLineCount:
  ClearErrors
   FileRead $UninstLog $R0
   IntOp $R1 $R1 + 1
   IfErrors 0 GetLineCount
 
 LoopRead:
  FileSeek $UninstLog 0 SET
  StrCpy $R2 0
  FindLine:
   FileRead $UninstLog $R0
   IntOp $R2 $R2 + 1
   StrCmp $R1 $R2 0 FindLine
 
   StrCpy $R0 $R0 -2
   IfFileExists "$R0\*.*" 0 +3
    RMDir $R0  #is dir
   Goto +3
   IfFileExists $R0 0 +2
   Delete $R0 #is file
 
  IntOp $R1 $R1 - 1
  StrCmp $R1 0 LoopDone
  Goto LoopRead
 LoopDone:
 FileClose $UninstLog
 Delete "$INSTDIR\${UninstLog}"

 Pop $R2
 Pop $R1
 Pop $R0
 
 Delete "$INSTDIR\Uninstall.exe"
 RMDir "$INSTDIR"
 
SectionEnd
