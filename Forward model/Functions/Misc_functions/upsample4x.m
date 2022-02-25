
function [I_upsampled4x] = upsample4x(I)
% fft function handles
F = @(x) fftshift(fft2(x));
Ft = @(x) ifft2(ifftshift(x));

% get image size
[m,n] = size(I);
% compute fft of image
I_fft = F(I);
% do zero padding
I_fft_padded = padarray(I_fft,[(3/2)*m (3/2)*n]);
% inverse fft of padded fft
tmp = Ft(I_fft_padded);
tmp = abs(tmp);
I_upsampled4x = double(tmp*(numel(tmp)/numel(I)));

end