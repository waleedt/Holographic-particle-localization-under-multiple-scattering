% Waleed Tahir
% 2017-04-04
% returns: 3D datacube for the simulation
function f = sim_obj(nx,ny,nz,n_obj,n_medium,num_particles)

    f = zeros(nx,ny,nz);
    f(nx/2,ny/2,1) = n_obj^2 - n_medium^2;
    f(nx/2,ny/2,nz) = n_obj^2 - n_medium^2;

    % create num_particles # of bubbles in a 256x256x20 cube
    fprintf('total no. of particles: %d\n',num_particles);
    %f_cent = zeros(nx/4,ny/4,nz);
    %bubble_idx = randi((nx/4)*(ny/4)*nz,[1 num_particles]);
    %f_cent(bubble_idx) = (n_obj^2)-n_medium^2;

    % put the 256x256x20 cube above in the middle of the 1024x1024x20 cube
    f = zeros(nx,ny,nz);
    %f(384:639,384:639,:) = f_cent;
    
    bubble_idx = randi(nx*ny*nz,[1 num_particles]);
    f(bubble_idx) = (n_obj^2)-n_medium^2;
end