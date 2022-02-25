% Waleed Tahir
% 2017-04-04
% Generate parameters to be used in the simulation

% parameters: wavelength and refractive indices
n_medium = 1.33;
n_obj = 1;
lambda=0.658/n_medium; % (um)

% parameters: da*nzta size
nx = 500;
ny = nx;
nz = 30;
deltaZ = 34.48; % um
offsetZ = 10*1000; % um (12mm)

dpix = 10;  % um
Xlen = 500*dpix;
Ylen = 500*dpix; % (883.2 um)


deltaX = Xlen/nx; % um
deltaY = Ylen/ny; % um


% this was a test for something, i cant remember what
less_than_pi = (pi*2*deltaX*(Xlen/2))/(lambda*offsetZ);




