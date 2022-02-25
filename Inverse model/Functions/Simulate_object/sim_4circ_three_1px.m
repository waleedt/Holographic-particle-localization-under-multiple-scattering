% return object comprising of four slices as mentioned below

function f = sim_4circ_three_1px(nx,ny,nz,n_obj,n_medium,~)


bndry_offset_px = 50;
circ_size_px = 6;
slice_size_px = nx;
numcircles = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% slice 1
% generate random indexes for circles
idx_range = [bndry_offset_px,slice_size_px-bndry_offset_px];
idx_x = randi(idx_range,[1 numcircles]);
idx_y = randi(idx_range,[1 numcircles]);
r = (circ_size_px/2);
idx_x = [nx/2-r nx/2-r nx/2+r+1 nx/2+r+1]; %center
idx_y = [ny/2-r ny/2+r+1 ny/2-r ny/2+r+1]; %center
% remove indexes near boundary
idx_x(idx_x<2*circ_size_px) = 0;
idx_y(idx_y<2*circ_size_px) = 0;
% make indexes negative to positive instead of one to positive
idx_x = idx_x - slice_size_px/2;
idx_y = idx_y - slice_size_px/2;

size = slice_size_px;
diameter = circ_size_px/size;

value = ones(numcircles,1);
vertlen = diameter * ones(numcircles,1);
horzlen = diameter * ones(numcircles,1);
xcor = idx_x';
ycor = idx_y';
phi = zeros(numcircles,1);

E = [value vertlen horzlen (xcor/size)*2 (ycor/size)*2 phi];

f = zeros(nx,ny,nz);
f(:,:,1) = phantom(E, size);
f(:,:,3) = phantom(E, size);



%% refractive index value
f = f * (n_obj^2 - n_medium^2);

%% debug
debug = 1;
if(debug)
   figure(1)
   for i = 1:nz
      subplot(1,nz,i)
      imagesc(f(:,:,i));
      colorbar;
      axis image;
   end
end

end
