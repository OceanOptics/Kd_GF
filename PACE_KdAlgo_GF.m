% Kd algorithm for PACE 

% Requires the following inputs that are listed below: 
%[TOTSCATAU,TOTEXTTAU,PS,a_aer,TOTANGSTR,wv_ang,lambda,sat_dt,lat,lon,a,bb]

% ---Ancillary data ----
% Scattering Efficiency at 550 nm : TOTSCATAU; 
% Extinction Efficiency at 550 nm : TOTEXTTAU;
% Atmospheric Pressure : PS;
% Anisotropy of phase function (Constant): a_aer;
a_aer = 0.67;
% Inputs for Tau mol calculation(Constant)
A = 85.35 *10^-4; B = -1.225 * 10^-4; C = 1.40 * 10^-4;

% --Level 2 data OR Ancillary Data ---
% Angstrom coefficient : TOTANGSTR;
% Reference wavelength for the Angstrom coefficient: wv_ang;
% wavelength of interest: lambda;
% Datetime of satellite overpass : sat_dt;
% Latitude of pixel : lat;
% Longitude of pixel: lon; 
% absorption(lambda) - Retrieved from QAA or any inversion algorithm : a;
% backscattering(lambda) from inversion algorithm : bb;


% --------- Computation of inputs for f  -------------------------------

%  Solar zenith angle in the air and in water (In radians)
theta_air = (sun_position_GB(sat_dt, lat, lon,0,0)) .* pi./180; % Angle at time of satellite overpass
theta_water = asin(sin(theta_air)./1.34); % Angle in the water

%  Molecular optical thickness : Tau mol
    % Requires PS in Mbar and lambda in microns - the following is converting
    % the Merra output (in Pascal) and wavelength in nanometers
tau_mol = (A./(lambda*10^-3)^4 )+ B./(lambda*10^-3)^5 +C/(lambda*10^-3)^6 *...
              PS*10^-2./1013.25;

%  Aerosols optical thickness : Tau aer 
tau_aer = TOTANGSTR .* (wv_ang/lambda).^ -TOTANGSTR;
          
%  Single Scattering Albedo : SSA_lambda
SSA_550 = TOTSCATAU / TOTEXTTAU;
SSA_lambda =  SSA_550 .* (550./lambda).^ -TOTANGSTR;

% ----------- Computation of f ------------------------------------------

%  Direct Transmittance (T_dir) 
t_dir = exp(-(tau_mol + tau_aer) ./ cos(theta_air));


%  Total Transmittance (T_tot)
t_tot = exp(-tau_mol .*0.5 ./ cos(theta_air)) .* ...
    exp( - (( 1- SSA_lambda .* 0.5 .* (1 + a_aer) ) .* tau_aer) ./ cos(theta_air));

% f 
f_lambda = t_dir./t_tot; 

% ------- Computation of Kd(lambda)-------------------------------------

D0_lambda = f_lambda./cos(theta_water) + 1.197 *(1 - f_lambda);

Kd_lambda = (a_lambda + bb_lambda).* D0_lambda; 


