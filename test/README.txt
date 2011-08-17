
% Configuration Files

 param.txt          - Sample Arepo parameter file.
 Config_ArepoVTK.sh - Sample Arepo configuration file.

% ArepoVTK Debugging Float Format Tests
 
 test.txt  - 2^3 constant rho=1.1 field [0,0,0]-[1,1,1].
 test1.txt - 101^3 constant rho=1.1 field, with rectangular prism overdensity 
             and spherical underdensity [0,0,0]-[1,1,1].

% Arepo HDF5 IC Format Tests

 Arepo1  - 2x2x1 (in z) uniform density grid, with one point at 2x density
 Arepo2  - 2^3 uniform rho grid, with one point at 2x density
 Arepo3  - 4^3 uniform rho grid, with one point at 2x density
 Arepo2b - 2^3 uniform rho grid with added central point (0.5,0.5,0.5) at 2x density

 *.hdf5.dump.txt - Output of ArepoMesh::DumpMesh() for the smaller ICs.

