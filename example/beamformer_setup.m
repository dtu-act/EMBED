% Get beamformer map
clear all
close all

load('micarray.mat');   % Load microphone positions
rn = pizza_array;
N = 50;    % Number of grid points in each dim
z0 = 5;     % Source distance [m]
phi = 15;   % Off-axis angle [degrees]
f = 1500;   % Frequency [Hz]
SNR = 15;   % Signal-to-noise ratio [dB]

source = int64([N/2-N/4 N/2; N/2+N/4 N/2]);    % x,y position of sources

[b,PSF,X,Y,x,y] = psf(N,z0,f,phi,rn,source,SNR);

figure
imagesc(real(b)), hold on
plot(source(:,1),source(:,2),'w*')
colormap(flipud(winter))

title('Normalized DAS beamforming map')
axis equal
colorbar