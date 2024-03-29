% Driver for miecoated.m, originally written by Christian matzler (see
% Matzler, 2002).
% inputs are:
%   rice : radius of inner ice sphere
%   rwater : radius of outer water sphere (from centre of ice sphere to
%           edge of water sphere)
%   
% Outputs: 
%   miecoated.m returns the following efficiency factors: [qext qsca qabs qb asy qratio]
%   the driver convolves these with the particle dimensions to return the
%   cross sections for extinction, scattering and absorption plus the
%   asymmetry parameter, q ratio and single scattering albedo.
%
% Note that the call to miecoated.m includes a term 'opt', which can be set
% to the values 1,2 or 3. These are alternative methods for calculating the
% absorption efficiency which can perform differently in some rare cases.
% The default is set to 3 here as this seems to perform better for large
% size parameters, which are likely to be used when modelling melting ice.
%
% Also note that the code includes an interpolation regime. This is because
% the original code produced NaNs for a few wavelengths at certain size
% parameters, particularly in the mid NIR wavelengths.


% Written by Joseph Cook, Feb 2017, University of Sheffield, UK

% set up variables
extinction=0;
scattering = 0;
absorption = 0;
asymmetry = 0;
q_ratio = 0;
ssa = 0;

% Define relevant parameters: Only rice and rwater are user defined.
rice = 1000; % inner sphere diameter in um
rwater = 1500; % outer sphere diameter in um (i.e. total coated sphere, not water layer thickness)

% calculate volume and density of sphere

XSArea_inner = pi*((rice)^2); % cross sectional area of ice core
XSArea_outer = pi*((rwater)^2) - XSArea_inner; % cross sectional area of water layer sphere
TotalXS = pi*((rwater)^2);


WatDensity = 999; % density of water at 1 degree C in kg m-3
IceDensity = 934; % density of ice in kg m-3

IceVol = 4/3 * pi * (rice)^3;
WatVol = (4/3 * pi * (rwater)^3) - IceVol;
TotalVol = 4/3 * pi * (rwater)^3;


IceMass = IceVol * IceDensity;
WatMass = WatVol * WatDensity;
TotalMass = IceMass + WatMass;


%define wavelength range
WL= 0.305:0.01:5;


for i = (1:1:470)
    % size parameters for inner and outer spheres
   x = 2 * pi * rice / (WL(i));

   y = 2 * pi * rwater / (WL(i));
   
    % read in refractive indices from refracICE and refracWATER per wavelength
   [RN CN] = refracICE(WL(i));
   m1 = complex(RN,CN);
   [RN2 CN2] = refracWATER(WL(i));
   m2 = complex(RN2,CN2);
   % call Miecoated and return efficiencies
   [qext qsca qabs qb asy qratio] = Miecoated(m1,m2,x,y,3);
   %append efficiencies to lists
   extinction(i)=qext;
   scattering(i) = qsca;
   absorption(i) = qabs;
   backscattering(i) = qb;
   asymmetry(i) = asy;
   q_ratio(i) = qratio;
   ssa(i) = qsca/qext;
end

% replace any possible nans with values estimated by cubic interpolation
extinction = fillmissing(extinction,'spline')
scattering= fillmissing(scattering,'spline')
absorption= fillmissing(absorption,'spline')
backscattering= fillmissing(backscattering,'spline')
asymmetry= fillmissing(asymmetry,'spline')
q_ratio= fillmissing(q_ratio,'spline')
ssa= fillmissing(ssa,'spline')

% calculate cross sections from efficiency factors
ExtXC = (extinction .* TotalXS);
ScaXC = (scattering .* TotalXS);
AbsXC = (absorption .* TotalXS);

ExtXCvol = (extinction.*TotalVol);
ScaXCvol = (scattering.*TotalVol);
AbsXCvol = (absorption.*TotalVol);

ExtXCmass = (extinction.*TotalMass);
ScaXCmass = (scattering.*TotalMass);
AbsXCmass = (absorption.*TotalMass);

% print fraction of sphere made up of water by mass and volume
water_frac_mss = 100-(IceMass/TotalMass)*100
water_frac_vol = 100-(IceVol/TotalVol)*100
ice_frac_mss = 100 - water_frac_mss
ice_frac_vol = 100 - water_frac_vol


part_dens = ((IceDensity*ice_frac_mss/100)+(WatDensity*water_frac_mss/100)) % density of particle as average weighted by mass of components

%plot relevant variables
% plot(WL,ExtXCmass)
% xlabel('Wavelength')
% ylabel('Cext')

% figure
plot(WL,ssa)
xlabel ('wavelength')
ylabel('ssa')
hold on
% 
% figure
% plot(WL,extinction)
% xlabel ('wavelength')
% ylabel('extinction efficiency')
% 
% figure
% plot(WL,scattering)
% xlabel ('wavelength')
% ylabel('scattering efficiency')
% 
% figure
% plot(WL,absorption)
% xlabel ('wavelength')
% ylabel('absorption efficiency')
% 
% figure
% plot(WL,asymmetry)
% xlabel ('wavelength')
% ylabel('asymmetry parameter (g)')