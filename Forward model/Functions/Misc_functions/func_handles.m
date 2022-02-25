% Defining FFT function handles
% F = @(x) fftshift(fft2(x));
% Fpad = @(x) fftshift(fft2(x,2*size(x,1),2*size(x,2)));
% Ft = @(x) ifft2(ifftshift(x));

F = @(x) fft2(x);
FG = @(x) fft2(ifftshift(x));
Fpad = @(x) fft2(x,2*size(x,1),2*size(x,2));
Ft = @(x) ifft2(x);
subindex = @(A, r, c) A(r, c);

F3D = @(x) fftn(x);
FG3D = @(x) fftn(ifftshift(x));
Fpad3D = @(x) fftn(x,[2*size(x,1),2*size(x,2),2*size(x,3)]);
Ft3D = @(x) ifftn(x);