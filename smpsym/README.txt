Tested on MSVC 2017 only (developed having in mind cross-platform use)


Requirements

* cinder (in windows: clone git repo, use cmake, build) + Cinder-ImGui
	https://github.com/cinder/Cinder

* libigl
	https://github.com/libigl/libigl


Optional

* geogram (for generation of triangulated stencil from input shape geometry)
    GEOGRAM_LIB_ONLY:BOOL=ON
    VORPALINE_PLATFORM:STRING=Darwin-clang
	[or else for other platforms, see geogram/cmake/platforms]

* integration with Matlab (for plotting)
	windows: add system PATH value matlabroot\extern\bin\win64
	see https://nl.mathworks.com/help/matlab/matlab_external/build-c-engine-programs.html


Other dependencies (included in folder external)

* cereal
	https://uscilab.github.io/cereal/

* eigen
	http://eigen.tuxfamily.org/index.php?title=Main_Page
	

Build

	Before calling cmake, copy file <CMakeLocal-template.txt> as <CMakeLocal.txt> and replace paths
	to the listed libraries in your system (cereal, cinder, libigl). Build the executable.


Usage

	1. Do one of the following:
		
		a. Drag and drop an OBJ file (stencil)
		
		OR
		
		b. Select in main menu: File -> Load simulation ...
			choose a BIN file from folder <shells>
	
	2. Start simulation


---------------------------------------------------------------------------------------------------

Possible issues with Cinder:
	- always build Cinder using cmake
	
Possible issues with Cereal
	- crashes on os x
		https://github.com/USCiLab/cereal/issues/439
