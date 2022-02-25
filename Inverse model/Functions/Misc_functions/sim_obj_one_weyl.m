% Waleed Tahir
% 2017-04-04
% returns: 3D datacube for the simulation
function f = sim_obj_one_weyl(nx,ny,nz,n_obj,n_medium)

    f = zeros(nx,ny,nz);
    f(nx/2,ny/2,1) = 1;

end