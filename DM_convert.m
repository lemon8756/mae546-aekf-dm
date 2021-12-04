% DM_convert.m

% Script to convert a DM voltage map to a DM surface map, and vice versa.

% Created on May 5 2015 by A.J. Eldorado Riggs

close all; clear all; clc;

D = 10e-3; % meters
Nact = 8; % Number of actuators across the DM
Np = 64/2; % Half the number of points across the pupil.
dx = D/2/Np;
xs = (-Np+1/2:1:Np-1/2)'*dx;
actFWHM = (D/Nact*1.4); %FWHM of an actuator, in meters

% Create the conversion matrix
infdx = zeros(Nact,length(xs));
for q = 1:Nact
    x_cent = q-Nact/2-1/2;
    infdx(q,:) = exp(-4*log(2)*((xs-D*x_cent/(Nact-1)).^2)./(actFWHM)^2);
end
infdx = 0.76*infdx;  % To give a peak value of 1 when all ones are input

figure; imagesc(infdx); axis equal xy tight;  title('Conversion Matrix'); colorbar;

%% Example: Go from a voltage map (u1) to a DM surface (h1)
u1 = magic(Nact); % ones(Nact);
h1 = infdx.'*u1*infdx;
figure; imagesc(u1);axis equal xy tight; title('Voltage Map 1'); colorbar;
figure; imagesc(h1);axis equal xy tight; title('Surface Map 1'); colorbar;

%% Example: Go from a DM surface h2 to a voltage map (u2)
close all;
% Generate a random phase map for h2:
% alpha=2*pi*rand(); [X,Y]=meshgrid(xs); k = 20;
% h2 = cos(2*pi*((Y+X*tan(alpha))/(sqrt(1+tan(alpha)^2)))*k+2*pi*rand()); 
h2 = kolmogorov(1,64,1/0.1);

u2 = (infdx*infdx')\infdx*h2*infdx'/(infdx*infdx'); % Fun times

figure; imagesc(h2);axis equal xy tight; title('Surface Map 2'); colorbar;
figure; imagesc(u2);axis equal xy tight; title('Voltage Map 2'); colorbar;
