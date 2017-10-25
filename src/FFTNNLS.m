function [x,info] = FFTNNLS(fun, PSF, b, x0, maxit)
% FFT-NNLS GRADIENT PROJECTION ALGORITHM WITH EXACT LINE SEARCH
%
% Usage:  [x,info] = FFTNNLS(@fun, PSF, b, x0, opt)
%
% Input:
%   @fun: Function handle with objective function and gradient
%   PSF: Point-spread function
%   b : Beamformer map
%   x0: Starting vector
%   maxit: Maximum number of iterations
%
% Output: 
%   x: Source distribution for all iterations
%   info: Struct with info about
%       .obj: Objective function values as a function of iterations
%       .time: Total time of algorithm
%
% Author: Oliver Lylloff
% Date: 25/9/14
% Latest revision: 25/9/14
%
% 
% Reference: 
% Ehrenfried, Klaus, and Lars Koop. 2008. 
% 'A Comparison of Iterative Deconvolution Algorithms for the Mapping of Acoustic Sources.' 
% American Institute of Aeronautics and Astronautics.

start_time = tic;

x = x0;

% Precompute fft of PSF
Fps = fft2(PSF);
FpsT = fft2(rot90(PSF,2));

fgx = @(x) fun(PSF,b,x,Fps,FpsT);      
n = 0;
while n < maxit
    n = n+1;
    [f,grad,r] = fgx(x);        
    d = grad;                   
    d(x == 0 & d > 0) = 0;      
    
    g = fftshift(ifft2(fft2(d).*Fps));         
    t = dot(g(:),r(:))/dot(g(:),g(:));    
    
    x = max(0,x - t*d);
    
    info.obj(n) = f;
end
info.time = toc(start_time);
end