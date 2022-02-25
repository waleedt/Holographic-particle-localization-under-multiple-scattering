% pixel sampling function
function[px_sampled_holo] =  pixel_sampling(holo,super_res_factor)

% filter kernel
h = (1/super_res_factor^2)*ones(super_res_factor,super_res_factor);
% filtering via 2D convoltuion
holo_filtered = conv2(holo,h,'valid');
% downsampling
px_sampled_holo = holo_filtered(1:super_res_factor:end,1:super_res_factor:end);

end