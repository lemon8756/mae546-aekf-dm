% kolmogorov.m

% Created on May 4, 2015 by Aaron James Lemmer

function [phsout,RMSphs,x] = kolmogorov(Dap,nxy,r0)

%Generate spectrum
dxy = Dap/nxy; 
df = 1/(2*nxy*dxy);
x = (-(nxy)/2:(nxy)/2-1)*dxy;
xf = (-(nxy)/2:(nxy)/2-1)*df;
[XF,YF]=meshgrid(xf,xf);
RF = sqrt(XF.^2+YF.^2);

%Calculate spectral density
SD = (0.0229/r0^(5/3)) .* (RF).^(-11/3);
SD(nxy/2+1,nxy/2+1)=0; 

% figure(1);
% imagesc(xf,xf,log10(SD),[-10,0]);
% axis square;
% colorbar;
% axis xy;
% xlabel('Spatial Frequency [1/m]');
% ylabel('Spatial Frequency [1/m]');

%Generate random phase screen
rndScr = randn(nxy,nxy) + 1j.*randn(nxy,nxy);
phs = real(ifft2(fftshift(sqrt(SD) .* rndScr)*nxy*nxy));

phsout = phs - mean(phs(:));
% figure(2);
% imagesc(x,x,(phs),[-1,1]);
% xlabel('[m]');
% ylabel('[m]');
% colorbar;
% axis square;
% axis xy;
RMSphs = sqrt(mean(phs(:).^2));
% RMS_notilt_ideal=sqrt(.134) * ((Dap/r0).^(5/6))
% RMS_uncomp_ideal=sqrt(1.02) * ((Dap/r0).^(5/6))