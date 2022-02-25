%% save uk
scatt_uk_dir = [obj_holo_dir,'\\scatt_field\\uk_total'];
mkdir(scatt_uk_dir);

%% save mat file of uk
uk_scatt_matfile = [scatt_uk_dir,'\uk_fields_all'];
save(uk_scatt_matfile,'E_rb_all');

%% plot and save some statistics of the norm of uk

% l2 norm of 3D field
uk_norm = zeros(length(born_orders));
for l = 1:length(born_orders)
    uk = E_rb_all(:,:,:,born_orders(l));
    uk_norm(l) = norm(uk(:), 2);
end
gca7 = figure(7);
plot(born_orders,uk_norm,'-*');
for l=1:length(born_orders)
    strpts = num2str(uk_norm(l));
    text(born_orders(l),uk_norm(l),strpts);
end
title('L2 norm of multiple scattered 3D field in the object vs Born order');
xlabel('Born order');
ylabel('Performance metric');
grid on;
filename =[scatt_uk_dir,'\3Dfields_l2norm_',num2str(num_particles(i)),'.png'];
saveas(gca7,filename,'png');

% l2 norm of 3D field slice-wise
uk_norm_2D = zeros(length(born_orders),nz);
for l = 1:length(born_orders)
    for ll = 1:nz
        uk_2D = E_rb_all(:,:,ll,l);
        uk_norm_2D(l,ll) = norm(uk_2D(:), 2);
    end
end
gca9 = figure(9);
plot(1:nz,uk_norm_2D(1,:),'-*');
hold on;
for l = 2:length(born_orders)
    plot(1:nz,uk_norm_2D(l,:),'-*');
end
grid on;
hold off;
xlabel('slice number');
ylabel('l2 norm of 2D slice');
title('l2 norm of 2D slices of uk, vs slice number, for all born orders');
cells=num2cell(born_orders);
for l = 1:length(cells)
   cellstr{l} = ['born',num2str(cells{l})]; 
end
legend(cellstr);

filename =[scatt_uk_dir,'\2Dfields_l2norm_',num2str(num_particles(i)),'.png'];
saveas(gca9,filename,'png');

%% save images of uk
scatt_uk_re_dir = [scatt_uk_dir,'\\real'];
scatt_uk_im_dir = [scatt_uk_dir,'\\imag'];
scatt_uk_abs_dir = [scatt_uk_dir,'\\abs'];
mkdir(scatt_uk_re_dir);
mkdir(scatt_uk_im_dir);
mkdir(scatt_uk_abs_dir);

%% save z-slice-wise
for l = 1:length(born_orders)
    scatt_uk_re_zw_dir{l} = sprintf('%s\\z-slice-wise\\born_%d',scatt_uk_re_dir,born_orders(l));
    mkdir(scatt_uk_re_zw_dir{l});
    scatt_uk_im_zw_dir{l} = sprintf('%s\\z-slice-wise\\born_%d',scatt_uk_im_dir,born_orders(l));
    mkdir(scatt_uk_im_zw_dir{l});
    scatt_uk_abs_zw_dir{l} = sprintf('%s\\z-slice-wise\\born_%d',scatt_uk_abs_dir,born_orders(l));
    mkdir(scatt_uk_abs_zw_dir{l});
end

gca8 = figure(8);
set(gca8,'position',[0 0 1224 1024])
% real
for l = 1:length(born_orders)
    real_cube = real(E_rb_all(:,:,:,l));
    maxval = max(real_cube(:));
    minval = min(real_cube(:));
    for ll = 1:nz
        imagesc(real(E_rb_all(:,:,ll,l)));
        colorbar
        if (minval ~= maxval)
            caxis([minval maxval]);
        end
        title(sprintf('real scatt field,born:%d,slice:%d',born_orders(l),ll));
        filename =[scatt_uk_re_zw_dir{l},'\u',num2str(born_orders(l)),'_re_slice_',num2str(ll),'.png'];
        saveas(gca8,filename,'png');
    end
end
% imag
for l = 1:length(born_orders)
    imag_cube = imag(E_rb_all(:,:,:,l));
    maxval = max(imag_cube(:));
    minval = min(imag_cube(:));
    for ll = 1:nz
        imagesc(imag(E_rb_all(:,:,ll,l)));
        colorbar
        if (minval ~= maxval)
            caxis([minval maxval]);
        end
        title(sprintf('imag scatt field,born:%d,slice:%d',born_orders(l),ll));
        filename =[scatt_uk_im_zw_dir{l},'\u',num2str(born_orders(l)),'_im_slice_',num2str(ll),'.png'];
        saveas(gca8,filename,'png');
    end
end
% abs
for l = 1:length(born_orders)
    abs_cube = abs(E_rb_all(:,:,:,l));
    maxval = max(abs_cube(:));
    minval = min(abs_cube(:));
    for ll = 1:nz
        imagesc(abs(E_rb_all(:,:,ll,l)));
        colorbar
        if (minval ~= maxval)
            caxis([minval maxval]);
        end
        title(sprintf('abs scatt field,born:%d,slice:%d',born_orders(l),ll));
        filename =[scatt_uk_abs_zw_dir{l},'\u',num2str(born_orders(l)),'_abs_slice_',num2str(ll),'.png'];
        saveas(gca8,filename,'png');
    end
end

%% save born-order-wise
for ll = 1:nz
    scatt_uk_re_bw_dir{ll} = sprintf('%s\\born-order-wise\\slice_%d',scatt_uk_re_dir,ll);
    mkdir(scatt_uk_re_bw_dir{ll});
    scatt_uk_im_bw_dir{ll} = sprintf('%s\\born-order-wise\\slice_%d',scatt_uk_im_dir,ll);
    mkdir(scatt_uk_im_bw_dir{ll});
    scatt_uk_abs_bw_dir{ll} = sprintf('%s\\born-order-wise\\slice_%d',scatt_uk_abs_dir,ll);
    mkdir(scatt_uk_abs_bw_dir{ll});
end

%gca8 = figure(8);
% real
for ll = 1:nz
    real_cube = real(E_rb_all(:,:,ll,:));
    maxval = max(real_cube(:));
    minval = min(real_cube(:));
    for l = 1:length(born_orders)
        imagesc(real(E_rb_all(:,:,ll,l)));
        colorbar
        if (minval ~= maxval)
            caxis([minval maxval]);
        end
        title(sprintf('real scatt field,born:%d,slice:%d',born_orders(l),ll));
        filename =[scatt_uk_re_bw_dir{ll},'\u',num2str(born_orders(l)),'_re_slice_',num2str(ll),'.png'];
        saveas(gca8,filename,'png');
    end
end
% imag
for ll = 1:nz
    imag_cube = imag(E_rb_all(:,:,ll,:));
    maxval = max(imag_cube(:));
    minval = min(imag_cube(:));
    for l = 1:length(born_orders)
        imagesc(imag(E_rb_all(:,:,ll,l)));
        colorbar
        if (minval ~= maxval)
            caxis([minval maxval]);
        end
        title(sprintf('imag scatt field,born:%d,slice:%d',born_orders(l),ll));
        filename =[scatt_uk_im_bw_dir{ll},'\u',num2str(born_orders(l)),'_im_slice_',num2str(ll),'.png'];
        saveas(gca8,filename,'png');
    end
end
% abs
for ll = 1:nz
    abs_cube = abs(E_rb_all(:,:,ll,:));
    maxval = max(abs_cube(:));
    minval = min(abs_cube(:));
    for l = 1:length(born_orders)
        imagesc(abs(E_rb_all(:,:,ll,l)));
        colorbar
        if (minval ~= maxval)
            caxis([minval maxval]);
        end
        title(sprintf('abs scatt field,born:%d,slice:%d',born_orders(l),ll));
        filename =[scatt_uk_abs_bw_dir{ll},'\u',num2str(born_orders(l)),'_abs_slice_',num2str(ll),'.png'];
        saveas(gca8,filename,'png');
    end
end