EMBED :: Efficient Methods for BEamforming Deconvolution.

EMBED is a small MATLAB package accompanying the journal paper:

Oliver Lylloff, Efrén Fernández-Grande, Finn Agerkvist, Jørgen Hald, Elisabet Tiana Roig and Martin S. Andersen.
"Improving the efficiency of deconvolution algorithms for sound source localization."
J. Acoust. Soc. Am. 138, 172 (2015); http://dx.doi.org/10.1121/1.4922516

The package contains two folders:
src/
	- FFTNNLS.m : FFT-NNLS algorithm
	- FISTA.m : FISTA algorithm
	- nnls.m : nnls problem f(x) = 1/2*||Ax-b||_2^2
	- zeropad.m : zeropads arrays before deconvolution procedure
example/
	- beamformer_setup.m : setup of grid size, opening angle, SNR
	- psf.m : calculates point-spread function
	- deconv.m : setup of deconvolution algorithms with number of iterations and zeropadding
	- micarray.mat : coordinates (x,y) of microphone array
	- run_example.m : run example and plot results
	
Usage: Add src/ to path and run run_example.m
