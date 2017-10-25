function [x,info] = FISTA(fun, PSF, b, x0, maxit)
% FISTA FAST ITERATIVE SHRINKAGE-THRESHOLDING ALGORITHM
% 
% Usage:  [x,info] = FISTA(@fun, PSF, b, x0, maxit)
%
% Input:
%   @fun: Function handle with objective function and gradient
%   PSF: Point-spread function
%   b : Beamformer map
%   x0: Starting vector
%   maxit: Maximum number of iterations
%
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
% Beck, Amir, and Marc Teboulle. 2009. 
% 'A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems.' 
% SIAM Journal on Imaging Sciences 2 (1) (January): 183-202.
%
start_time = tic;

% Initialize variables
x = x0;
xold = x;
y = x;
t = 1; 

% Precompute fft of PSF
Fps = fft2(PSF);
FpsT = fft2(rot90(PSF,2));

% Evaluate f(x0)
fgx = @(x) fun(PSF,b,x,Fps,FpsT);      

% Compute Lipschitz constant
L = lipschitz(PSF,Fps);

% For n = 1
[fy,grady] = fgx(y);
info.obj(1) = fy;
% Start iteration
n = 0;
while n < maxit    
    
    n = n+1;
    x = max(0,y - (1/L)*grady);
    
    tnew = (1+sqrt(1+4*t*t))/2;
    y = x + ((t-1)/tnew)*(x-xold);
    
    [fy,grady] = fgx(y);
    xold = x;
    t = tnew;
    info.obj(n+1) = fy;     
end
info.time = toc(start_time);
end


function L = lipschitz(PSF,Fps)
% Estimate Lipschitz constant by power iteration
x = rand(size(PSF));
for k = 1:10
    x = fftshift(ifft2(fft2(x).*Fps))/norm(x,'fro');
end
    L = norm(x,'fro')^2;    % lambda(A'A) Assuming a symmetric matric A
end
