# Kd_GF
Compute the diffuse attenuation coefficient of light in water from Remote Sensing Reflectance (Kd) following Begouen et al., 2024. Hyperspectral, will work for any lambda 

Inputs are the following from Standard Level 2 Ocean Color Sensor: 
- lambda : wavelength of interest for which want Kd 
- a & bb retrieved from IOP inversion such as QAA.
- datetime of rrs measurement
- Latitude & longitude of rrs pixel
- wavelength of reference for angstrom coefficient 

Ancillary inputs needed (Can be downloaded from MERRA or in L1-L2 processing line) : 
Scattering Efficiency at 550 nm : TOTSCATAU; 
Extinction Efficiency at 550 nm : TOTEXTTAU;
Atmospheric Pressure : PS;

The code needsthe sun_position_GB.m function to compute solar angle in atmosphere.
