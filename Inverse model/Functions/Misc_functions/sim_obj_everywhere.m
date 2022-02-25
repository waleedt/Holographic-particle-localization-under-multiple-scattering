% Waleed Tahir
% 2017-04-04
% returns: 3D datacube for the simulation
function [f,bubble_idx] = sim_obj_everywhere(nx,ny,nz,n_obj,n_medium,lambda,num_particles)
    % initialize object
    f = zeros(nx,ny,nz);
    % randomly generate indices for num_particles number of bubbles
    bubble_idx = randi(nx*ny*nz,[1 num_particles]);
    % assign refractive indexes to bubbles
    k = 2*pi/lambda;
    f(bubble_idx) = (k^2/(4*pi))*(n_obj^2 - n_medium^2);
end