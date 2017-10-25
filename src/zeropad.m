function Bpad = zeropad(B)
%ZEROPAD Pad array with zeros to avoid wrap-around effects in decovolution
%
%function Bpad = zeropad(B)
%
%      Bpad = zeropad(B);
%
%  Pad B (N by N) with zeros to make it 2N by 2N. 
%
%  Input:
%     B     Array to be zero-padded
%
%  Output:
%        Bpad  Padded 2N by 2N array.
%
% Reference: See Chapter 4,
%            "Deblurring Images - Matrices, Spectra, and Filtering"
%            by P. C. Hansen, J. G. Nagy, and D. P. O'Leary,
%            SIAM, Philadelphia, 2006.
%
% Oliver Lylloff
% Date: 25/9/2014
% Revision: 25/9/2014

[M,N] = size(B);
Bpad = zeros(M+N);
Bpad(int64(M/2)+1:int64(M/2 + M),int64(M/2)+1:int64(M/2 + M)) = B;
end
