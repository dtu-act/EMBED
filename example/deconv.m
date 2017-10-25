% Run deconvolution algorithms
N = size(b,1);
b = real(zeropad(b));
PSF = zeropad(PSF);
x0 = zeros(2*N);
maxit = 3000;      % Set number of iterations

% FFT-NNLS:
[x_nnls,info_nnls] = FFTNNLS(@nnls,PSF,b,x0,maxit);
% Remove zero-padding
x_nnls = x_nnls(int64(N/2)+1:int64(N/2 + N),int64(N/2)+1:int64(N/2 + N));

% FISTA:
[x_fista,info_fista] = FISTA(@nnls,PSF,b,x0,maxit);
% Remove zero-padding
x_fista = x_fista(int64(N/2)+1:int64(N/2 + N),int64(N/2)+1:int64(N/2 + N));