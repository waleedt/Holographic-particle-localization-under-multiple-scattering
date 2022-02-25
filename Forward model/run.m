%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate datacube for holography
% Waleed Tahir
% Date: 3 Jan 2018
add_paths;

for px_spacing = 8:8
    
    n_obj = 1.43;
    for numspheres = [4]
        
        nx = 256;
        ny = 256;
        nz = 37;
        b_offset = 20;
        voxel_size = 3.45/20; %um
        
        particle_dia = 1/voxel_size;
        particle_rad = ceil(particle_dia/2);
        
        f = zeros(nx,ny,nz);
        tmp = zeros(nx,ny,nz);
        x = (1:nx) - nx/2;
        y = (1:ny) - ny/2;
        z = (1:nz) - (nz+1)/2;
        [X,Y,Z] = meshgrid(x,y,z);
        
        
        z_1slice(1) = z(particle_rad+1);
        z_1slice(2) = z(end-(particle_rad));
        
        % generate 3D object with random spheres
        for j = 1 : length(z_1slice)
            z_tmp = z_1slice(j);
            % generate indices of random sphere locations
            x_boundary = x(b_offset:end-b_offset);
            y_boundary = y(b_offset:end-b_offset);
            px = px_spacing;
            x_rand = [px px -px -px];
            y_rand = [px -px px -px];
            for i = 1:numspheres
                tmp(sqrt((X-x_rand(i)).^2+(Y-y_rand(i)).^2+((Z-z_tmp).^2)) <= particle_rad) = 1;
                f = f + tmp;
                f(f >= 1) = 1;
            end
        end
        
        res_folder = sprintf('../IM/Data/px_spacing_%d',px_spacing);
        for i = 1:7
            h = imagesc(f(:,:,i)); axis image;
            fn = sprintf('%s/n%1.2f/slices',res_folder,n_obj);
            mkdir(fn);
            filename = sprintf('%s/slice%d.png',fn,i);
            saveas(h,filename,'png');
        end
        
        fn = sprintf('%s/f_obj_1.mat',res_folder);
        f_obj = f;
        save(fn,'f_obj');
        
    end
    
    fm_born10_fn(n_obj,res_folder)
    
end