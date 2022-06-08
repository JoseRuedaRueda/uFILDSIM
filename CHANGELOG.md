## 2.1 Minor print
- Print also the number of markers which collide backward with the plate

## 2.0 Extra information in files
- Strike points files now include extra information on the applied models

## 1.6 Self-shadowing
- Included self-shadowing option in the FILDSIM simulation mode
- Include the filter to skip markers colliding with the back of the foil in the signal part of the code

## 1.5 3D bug fix?
- The 3D interpolation problem may have been solved

## 1.4 Minor upgrade:
- Time point is saved at the end of the signal file

## 1.3 StrikePoints:
- Added a flag to save the strike points in the scintillator. This allow to ignore that file if we just want the map, which saves time in the database scans

## 1.2 Bug solve:
- Solved a bug by which collimator strikes points could not be loaded in the Scintillator Suite
- Solved a bug in the 3D interpolation

## 1.1 Scintillator modelling:
- Models to include the yield in the scintillator added to the code
- FIDASIM weight now is scaled by the pinhole area so is stored in part/s
- Initial and final energy stored in the signal strike points file

## 1.0 CarbonFoil modeling:
- No needed to create the /bin/ forlder during the first installation
- Included different models of energy loss and attenuation in the carbon foil (for INPA)
- If the scintillator is not the last element of the geometry files, the code will internally re-arange it
- nGeomElements was moved to the ExtraGeometryParams namelist, as it is a characteristic of the chosen geometry
- Particle impinging from behind the scintillator will no longer be considered
- ScintNormal not needed in the Extra geometry params namelist, as it will be calculated from the triangles of the scintillator. Notice: The order of the points in the triangles should be right-handed such as the normal goes in the direction of the emmited light
- Number of geometry element to load is no longer needed to be be included in the configuration namelist, this parameter was moved to the Geometry namelist
- #bug#: maxT needed to be negative is backtracing mode was activated.

## 0.3 Namelist update:
- Default values are defined for namelist parameters, no need of giving all parameters like inputs
- GeomID changed to geomFolder: Now we would need to give the program the whole path to the folder with the element files and not just the geomID
- SINPA_dir changed to runFolder: Now instead the root directory of the SINPA code, the FORTRAN code need the folder containing the run (the parent folder for the inputs/results) This is included to gain flexibility with the SINPA execution to remap FILD frames
- Gyroradius and XI arrays are sorted in the FORTRAN code, such that there is no 'issues' with their reading in python

## 0.2 Backtracing:
- Included a flag to allow for backtracing in the namelist (the flag is called backtrace). The backtrace is just included reversing the sign of the time step in the integration
- Collimator factor definition now takes into account the restricted interval for launching the particles
