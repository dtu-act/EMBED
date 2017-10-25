% Results
clear all, close all
beamformer_setup;   % Get beamformer map
deconv;             % Run deconvolution algorithms

% Print time spent by algorithms
sprintf('Total time of algorithms after %d iterations:\n FFT-NNLS:\t%0.3g seconds \n FISTA:\t%0.3g seconds\n Speed-up: %0.3g X',maxit,info_nnls.time,info_fista.time,info_nnls.time/info_fista.time)

% Plot source distributions
figure
imagesc(x_nnls), hold on
colormap(flipud(winter))
plot(source(:,1),source(:,2),'w*')
title(['FFT-NNLS: Solution after ' num2str(maxit) ' iterations.'])
hold off

figure
imagesc(x_fista), hold on
colormap(flipud(winter))
plot(source(:,1),source(:,2),'w*')
title(['FISTA: Solution after ' num2str(maxit) ' iterations.'])

% Convergence 
niter = int64(maxit/3);
figure
loglog(1:niter,abs(info_nnls.obj(1:niter)-info_fista.obj(maxit))/info_fista.obj(maxit),'r'),hold on
loglog(1:niter,abs(info_fista.obj(1:niter)-info_fista.obj(maxit))/info_fista.obj(maxit),'b--')
xlabel('Iterations k','FontSize',14)
ylabel('f(x^k)-f(x^*)/f(x^*)','FontSize',14)