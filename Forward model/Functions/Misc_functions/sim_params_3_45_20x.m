% Waleed Tahir
% 2017-04-04
% Generate parameters to be used in the simulation

% parameters: wavelength and refractive indices
n_medium = 1.33;
lambda=0.6328/n_medium; % (um)

% parameters: da*nzta size
nx = nx;
ny = nx;
nz = rec_num_slices;
dpix = 3.45/20;  % um
if stype == 3 %sphere 
    deltaZ = dpix; % um (1.54mm)
else
    deltaZ = 8*(3.45/20); %um
end
offsetZ = 5; % um (17mm)

       
 
Xlen = nx*dpix; % (110.4 um)
Ylen = nx*dpix; % (110.4 um)


deltaX = Xlen/nx; % um
deltaY = Ylen/ny; % um

% NA of objective lens
NA = 0.4;