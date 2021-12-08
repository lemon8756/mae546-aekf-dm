function [phase,rmsPhase,x] = kolmogorov(effAper,nxy,r0)
%KOLMOGOROV   2-D Kolmogorov phase screen.
%   [phase,rmsPhase,x] = KOLMOGOROV(diamAper,nxy,r0) returns a zero-mean 
%   2-D scalar phase screen that models Kolmogorov turbulence, the RMS 
%   value of the phase screen, and the vector of coordinates that applies 
%   to both axes.
%
%   effAper: effective aperture, either the lateral width (square) or the
%            diameter
%       nxy: number of pixels in x and y directions
%        r0: Fried coherence length
%
%   KOLMOGOROV generates the phase screen as a convolution: the inverse
%   fast Fourier transform (FFT) of the product of the Wiener spectrum with
%   a filter of normally-distributed random complex numbers. This serves to
%   randomly select and scale components of the Wiener spectrum in the 
%   Fourier domain.
%
%   Although the Wiener spectrum is unbounded at the origin, the phase
%   difference between any two points across an aperture is bounded [1]. To
%   resolve this computationally, the coordinates are edge-defined such
%   that, for an even number of pixels, the origin lies at the intersection
%   of the boundaries of the center four pixels and there is no issue with 
%   the singularity.  For an odd number of pixels, the origin lies at the 
%   center of the middle pixel, so KOLMOGOROV sets this pixel to zero.
%
%   The FFT method utilized here is limited in that its sampling of the
%   spectrum at low spatial frequencies is too coarse and that the lowest
%   frequency components are neglected if the origin pixel is zeroed.
%
%   References:
%   [1] R. G. Lane, A. Glindemann, and J. C. Dainty. Simulation of a
%   Kolmogorov phase screen.  Waves in Random Media, 2(1992):209–224, 1992.
%   [2] D. L. Fried. Statistics of a geometric representation of wavefront 
%   distortion. J. Opt. Soc. Am., 55(11):1427–1435, 1965.

%   Created on May 4, 2015 by Aaron James Lemmer

dxy = effAper/nxy; 
df = 1/(2*nxy*dxy);
x = (-(nxy-1)/2:(nxy-1)/2)*dxy;
f = (-(nxy-1)/2:(nxy-1)/2)*df;
[fx,fy] = meshgrid(f,f);
spatialFreq = sqrt(fx.^2+fy.^2);

% Calculate spectral density using the Wiener spectrum
spectralDensity = (0.0229/r0^(5/3)) .* (spatialFreq).^(-11/3);
if mod(nxy,2) ~= 0
    % Set the center pixel to zero
    spectralDensity((nxy+1)/2,(nxy+1)/2) = 0;
end

% figure(1);
% imagesc(f,f,log10(spectralDensity),[-10,0]);
% axis square;
% colorbar;
% axis xy;
% xlabel('Spatial Frequency [1/m]');
% ylabel('Spatial Frequency [1/m]');

% Generate random filter to select and scale frequency components
rndFilter = randn(nxy,nxy) + 1j.*randn(nxy,nxy);

% Convolve the spectral density with the random filter
rndPhase = real(ifft2(fftshift(sqrt(spectralDensity).*rndFilter)*nxy*nxy));

phase = rndPhase - mean(rndPhase(:));
rmsPhase = sqrt(mean(phase(:).^2));

% figure(2);
% imagesc(x,x,(rndPhase),[-1,1]);
% xlabel('[m]');
% ylabel('[m]');
% colorbar;
% axis square;
% axis xy;