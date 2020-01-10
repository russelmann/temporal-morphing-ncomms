# Programming temporal morphing of self-actuated shells

Code and data for research paper Guseinov, R. et al. Programming temporal morphing of self-actuated shells. *Nature Communications* (2020).

### Abstract

Advances in shape-morphing materials, such as hydrogels, shape-memory polymers and light-responsive polymers have enabled prescribing self-directed deformations of initially flat geometries. However, most proposed solutions evolve towards a target geometry without considering time-dependent actuation paths. To achieve more complex geometries and avoid self-collisions, it is critical to encode a spatial and temporal shape evolution within the initially flat shell. Recent realizations of time-dependent morphing are limited to the actuation of few, discrete hinges and cannot form doubly curved surfaces. Here, we demonstrate a method for encoding temporal shape evolution in architected shells that assume complex shapes and doubly curved geometries. The shells are non-periodic tessellations of pre-stressed contractile unit cells that soften in water at rates prescribed locally by mesostructure geometry. The ensuing midplane contraction is coupled to the formation of encoded curvatures. We propose an inverse design tool based on a data-driven model for unit cells' temporal responses.

## Contents

* Folder **extra**: Processing and modeling for unit cells and shells. Plotting.
  - Folder **fabrication**: Helpers for fabrication
  - Folder **shell_mechanical_tests** : Shell mechanical tests data and code for plotting it
  - Folder **shell_scan**: Code for measuring shape precision via 3D scans
  - Folder **unit_cell_generator**: Generate unit cell specimens and tools for measurement setup
  - Folder **unit_cell_modeling**: Data and code for modeling unit cells, including plotting
  - File **plot_stress_relax.m**: Plot membrane stress relaxation
  - File **rubber_stress_relaxation.csv**: Membrane stress relaxation data
* Folder **shells**: Input shell geometries and simulation initialization data.
* Folder **smpsym**: Inverse design and simulation tool. (N.B.: for experienced C++ users only!)

## Inverse design and simulation tool (smpsym)

This code is tested only with MSVC 2017, but developed in a compatible way with gcc and Clang compilers.

### Requirements

* cinder (clone git repo, use cmake, build) + Cinder-ImGui
	https://github.com/cinder/Cinder

* libigl
	https://github.com/libigl/libigl

### Optional

* geogram (for generation of triangulated stencils)
    https://github.com/alicevision/geogram

* integration with Matlab (for plotting)

### Other dependencies (included in this repository)

* cereal
	https://uscilab.github.io/cereal/

* eigen
	http://eigen.tuxfamily.org/index.php?title=Main_Page
	
### Build

Before calling cmake, copy file <CMakeLocal-template.txt> as <CMakeLocal.txt> and replace paths
to the listed libraries in your system (cinder, libigl, etc.). Build the executable.

### Usage

To start editing *actuation time landscape* of a stencil represented by an OBJ file, drag and drop it to the app window.
Alternatively, to load a simulation file, select in main menu: **File** -> **Load simulation ...** and choose a BIN file, e.g. from folder **shells**.
Now you can edit *actuation time landscape* and start the simultaion.

## License
This code is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
[GNU General Public License](https://www.gnu.org/licenses/) for more details.
