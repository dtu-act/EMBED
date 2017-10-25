function [f,g,r] = nnls(PSF,b,x,varargin)
% NONNEGATIVE LEAST SQUARES OBJECTIVE FUNCTION
% 
% Compute objective, gradient and residual of 
% Nonnegative least squares problem (NNLS): 
%   minimize_x      0.5||Ax - b||^2 
%   subject to         x >= 0
%
% Usage:    
%       f = nnls(PSF,b,x)
%       [f,g] = nnls(PSF,b,x)
%       [f,g,r] = nnls(PSF,b,x)
%
% Input:
%   PSF: Point-spread function (matrix)
%   b : Beamformer map (matrix)
%   x: Source distribution (matrix)
%   varargin: Can take precomputed PSF
%
% Output: 
%   f: Objective function value (scalar)
%   g: Gradient (matrix)
%   r: Residual (matrix)
%
% Author: Oliver Lylloff
% Date: 25/9/14
% Latest revision: 25/9/14
%

if nargin==5
    Fps = varargin{1};
    FpsT = varargin{2};
else
    Fps = fft2(PSF);
    FpsT = fft2(rot90(PSF,2));
end

r = fftshift(ifft2(fft2(x).*Fps)) - b;
f = 0.5*norm(r,'fro')^2;
if (nargout > 1)
    g = fftshift(ifft2(fft2(r).*FpsT));
end

end