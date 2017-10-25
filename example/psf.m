function [b,PSF,X,Y,x,y] = psf(N,z0,f,phi,rn,source,SNR)

% Number of microphones in array
M = size(rn,1); 

% Constants
c = 343;               % Speed of sound
omega = 2*pi*f;
Np = round(N*1.3);     % PSF grid size

% Allocate memory
dj = zeros(Np,Np,M);
ej = zeros(Np,Np,M);
ejhat = zeros(Np,Np,M);
PSF = zeros(Np,Np);
b = zeros(N,N);

% Setup:
L = 2*z0*tand(phi);               % Area covered at distance z0

% Region of interest
x = [-L/2 L/2];     % x grid 
y = x;              % y grid

% Beamformer grid
rx = linspace(x(1),x(2),N);
[X,Y] = meshgrid(rx);

% PSF grid
rx_psf = linspace(x(1),x(2),Np);
[Xp,Yp] = meshgrid(rx_psf);

% Following [Koop et al.]: 
d0 = sqrt(Xp.^2 + Yp.^2 + z0^2);    % |d0|: Distance from (0,0,0) to all mesh points

% Distance dj from each microphone to each grid point 
% Steering vectors ej, ejhat
for n = 1:M
    dj(:,:,n) = sqrt((Xp-rn(n,1)).^2+(Yp-rn(n,2)).^2 + z0^2);
    ej(:,:,n) = (dj(:,:,n)./d0).*exp(1j*omega.*dj(:,:,n)./c);
    ejhat(:,:,n) = (d0./dj(:,:,n)).*exp(1j*omega.*dj(:,:,n)./c);
end

% POINT-SPREAD FUNCTION
% Is constructed by the DAS beamformer response of a unit strength point
% source in the middel of the grid
ind = round(Np/2);

e_unit = squeeze(ejhat(ind,ind,:));

for ii = 1:length(Xp)
    for jj = 1:length(Yp)
        PSF(ii,jj) = dot(squeeze(ej(ii,jj,:)),e_unit);
    end
end

PSF = rot90(abs(PSF).^2/M^2);      % Normalize PSF
if N ~= Np
    PSF = interp2(Xp,Yp,PSF,X,Y);
end

% CROSS-SPECTRAL MATRIX
dj = zeros(N,N,M);
ej = zeros(N,N,M);
ejhat = zeros(N,N,M);
% Distance and steering vectors
d0 = sqrt(X.^2 + Y.^2 + z0^2);    % |d0|: Distance from (0,0,0) to all mesh points
for n = 1:M
    dj(:,:,n) = sqrt((X-rn(n,1)).^2+(Y-rn(n,2)).^2 + z0^2);
    ej(:,:,n) = (dj(:,:,n)./d0).*exp(1j*omega.*dj(:,:,n)./c);
    ejhat(:,:,n) = (d0./dj(:,:,n)).*exp(1j*omega.*dj(:,:,n)./c)+10^(-SNR/10)*(rand(N,N)+1j*rand(N,N));
end

q = zeros(N,N);
Rmn = zeros(M,M);

for k = 1:size(source,2)
    q(source(k,2),source(k,1)) = 1;
    Rmn = Rmn + squeeze(ejhat(source(k,2),source(k,1),:))*squeeze(ejhat(source(k,2),source(k,1),:))';
end

% Diagonal removal
Rmn(logical(eye(size(Rmn)))) = 0;

% DELAY AND SUM BEAMFORMER OUTPUT
for ii = 1:length(X)
    for jj = 1:length(Y)
        e = squeeze(ej(ii,jj,:));
        b(ii,jj) = dot(e,Rmn*e)/(M^2-M);
    end
end
end
