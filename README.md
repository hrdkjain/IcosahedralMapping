This project contains matlab code for encoding manifold mesh to icosahedral grid using spherical parameterized representation of the mesh. This has been introduced in "GenIcoNet: Generative Icosahedral Mesh Convolutional Network", 3DV 2021.

### Usage
#### Download
```
git clone https://github.com/hrdkjain/Matlab-Functions.git
git clone https://github.com/hrdkjain/IcosahedralMapping.git
```
#### Environment Requirement 
The project has been tested on Ubuntu 20.04 with Matlab (version R2021a).
Using `orient_map.m` requires Simulink 3D Animation Toolbox for `vrrotvec2mat(...)` and Phased Array System Toolbox for `rotz(...)` amongst others.
#### Project Dependencies
This project depends on [Matlab-Functions](https://github.com/hrdkjain/m_libaries).
#### Execution
- Copy `include.m.example` to `include.m` and modify the paths. 
- List the manifold mesh and corresponding parameterized mesh using `fileLister.m`
- Run `main.m` to save encoded icosahedral mesh as .off and .mat files.

###  Projects
If you like this project you also might be interested in other project which use IcosahedralMapping for [mesh generation](https://github.com/hrdkjain/GenIcoNet).