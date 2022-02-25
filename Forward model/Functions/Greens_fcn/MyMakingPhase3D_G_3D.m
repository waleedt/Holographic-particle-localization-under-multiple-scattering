function Phase3D_G = MyMakingPhase3D_G_3D(Nx,Ny,Nz,lambda,...
    deltaX,deltaY,deltaZ,offsetZ,F,is_rm_alias,isconv_G,isMemOpt)

%% Zero padding for ani aliasing
if is_rm_alias
    pad_len = 2;
    Nx = Nx * pad_len;
    Ny = Ny * pad_len;
    Nz = Nz * pad_len;
end

%% Spatial frequency variables
k=1/lambda;
KX=[ceil(-Nx/2):1:ceil(Nx/2-1)]'.*(1/(Nx*deltaX));
KY=[ceil(-Ny/2):1:ceil(Ny/2-1)].*(1/(Ny*deltaY));
KZ=[ceil(-Nz/2):1:ceil(Nz/2-1)].*(1/(Nz*deltaZ));

[ky,kx,kz] = meshgrid(KX,KY,KZ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute via Conv approach
k=(2*pi)/lambda;
%% Spatial variables
% x-y-z axes coordinates
X_mat = kx * Nx * deltaX^2; % repeat X, 1 row Ny col
Y_mat = ky * Ny * deltaY^2; % repeat Y, Nx row 1 col
Z_mat = kz * Nz * deltaZ^2; % repeat Y, Nx row 1 col

%% Generate Phase3D_G
const = ((k^2)/(4*pi))*deltaX^3;

term = sqrt(X_mat.^2+Y_mat.^2+Z_mat.^2);

term(term==0)=0.1; % avoid Nan
Greens_spatial = exp(1i*k*term)./term;
Greens_spatial(term<0.172) = 0.0378 + 0.0544i;
Phase3D_G = F(Greens_spatial)*const;

disp('phase3D_G_conv done');

end

