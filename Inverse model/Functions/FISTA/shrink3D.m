function xhat = shrink3D(y, tau, isamp)
%%% Shrinkage function for complex 3D-TV denoising.
%%% Orig Author: U. S. Kamilov
%%% Edited by: Waleed Tahir
%%% Added capability to separately threshold real and imag parts of the
%%% signal. The thresholding method can be selected via a flag.

% shortcut for computing norm along the fourth dimension
computeNormInZ = @(x) sqrt(sum(abs(x).^2, 4));

% isamp : 1 -> apply thresholding to amplitude
% isamp : 0 -> apply thresholding separately to real amd imag parts
if isamp
    % convert data to a 6D array with real elements
    yv = zeros([size(y, 1), size(y, 2), size(y, 3), 6]);
    yv(:,:,:,1) = real(y(:,:,:,1));
    yv(:,:,:,2) = imag(y(:,:,:,1));
    yv(:,:,:,3) = real(y(:,:,:,2));
    yv(:,:,:,4) = imag(y(:,:,:,2));
    yv(:,:,:,5) = real(y(:,:,:,3));
    yv(:,:,:,6) = imag(y(:,:,:,3));
    
    % compute the norm in dimension 4
    norm_yv = computeNormInZ(yv);
    
    amp = max(norm_yv-tau, 0);
    norm_yv(norm_yv<=0) = 1; % to avoid division by 0
    
    xhatv = repmat(amp./norm_yv, [1, 1, 1, 6]) .* yv;
    
    xhat = zeros(size(y));
    xhat(:,:,:,1) = xhatv(:,:,:,1) + 1j*xhatv(:,:,:,2);
    xhat(:,:,:,2) = xhatv(:,:,:,3) + 1j*xhatv(:,:,:,4);
    xhat(:,:,:,3) = xhatv(:,:,:,5) + 1j*xhatv(:,:,:,6);
else
    % convert data to two 3D arrays, each with real elements
    yvr = zeros([size(y, 1), size(y, 2), size(y, 3), 3]);
    yvi = zeros([size(y, 1), size(y, 2), size(y, 3), 3]);
    yvr(:,:,:,1) = real(y(:,:,:,1));
    yvi(:,:,:,1) = imag(y(:,:,:,1));
    yvr(:,:,:,2) = real(y(:,:,:,2));
    yvi(:,:,:,2) = imag(y(:,:,:,2));
    yvr(:,:,:,3) = real(y(:,:,:,3));
    yvi(:,:,:,3) = imag(y(:,:,:,3));
    
    % compute the norms in dimension 4
    norm_yvr = computeNormInZ(yvr);
    norm_yvi = computeNormInZ(yvi);
    
    ampr = max(norm_yvr-tau, 0);
    ampi = max(norm_yvi-tau, 0);
    norm_yvr(norm_yvr<=0) = 1; % to avoid division by 0
    norm_yvi(norm_yvi<=0) = 1; % to avoid division by 0
    
    xhatvr = repmat(ampr./norm_yvr, [1, 1, 1, 3]) .* yvr;
    xhatvi = repmat(ampi./norm_yvi, [1, 1, 1, 3]) .* yvi;
    
    xhat = zeros(size(y));
    xhat(:,:,:,1) = xhatvr(:,:,:,1) + 1j*xhatvi(:,:,:,1);
    xhat(:,:,:,2) = xhatvr(:,:,:,2) + 1j*xhatvi(:,:,:,2);
    xhat(:,:,:,3) = xhatvr(:,:,:,3) + 1j*xhatvi(:,:,:,3);
end

end