z_min= min(Holo(:));
z_max = max(Holo(:));
holo_img = mat2gray(Holo,[z_min,z_max]);
fn_h =['holo_est.bmp'];
imwrite(holo_img,fn_h,'bmp')