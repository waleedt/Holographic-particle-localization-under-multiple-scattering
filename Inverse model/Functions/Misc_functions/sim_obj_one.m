% Waleed Tahir
% 2017-04-04
% returns: 3D datacube for the simulation
function f = sim_obj_one(nx,ny,nz,n_obj,n_medium,lambda,~)
    k = 2*pi/lambda;
    f = zeros(nx,ny,nz);
    %f(nx/2,ny/2,1) = (k^2/(4*pi))*(n_obj^2 - n_medium^2);
    %f(1,1,1) = (k^2/(4*pi))*(n_obj^2 - n_medium^2);
    
    f(nx/2,ny/2,1) = (n_obj^2 - n_medium^2);

end