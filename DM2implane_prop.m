%DM2implane_prop.m

%Created on May 6, 2015 13:48 by Aaron James Lemmer

function [imout,xout,contrast] = DM2implane_prop(phsin,nxy)

% DM2implane propagates an input phase aberration to image plane

% phsin = aberrated phase at shaped pupil plane
% nxy = number of pixels in input plane

% imout = intensity at image plane
% xout = vector of coordinates in image plane

%% Variable initialization

D = 0.01;  % size of aperture [m]
lambda = 635*10^-9;  % wavelength [m]
f = 1.524;  % lens focal length (=60") [m]
u = 32*f*lambda/D;  % output plane size
nx = nxy;  % number of pixels in output plane

%% Calculate the electric field

Eab = exp(1i.*2*pi*phsin);
SP = MakeNewRippleMask('N11.dat',nx/2);
Ein = SP.*(1 + Eab);

%% Fourier transform Ein, get electric field at FP, Eout

[Eout,xout,~] = ft(Ein, f, D, u, lambda, nx, 1);

im = abs(Eout).^2./max(max(abs(Eout).^2));  %intensity of Eout
[Xout,Yout] = meshgrid(xout,xout);
mask = ~Circ(Xout,Yout,u/4).*Circ(Xout,Yout,20*f*lambda/D).*~Ang(Xout,Yout,50);
imout = im.*mask;


contrast = mean(imout(mask == 1));

function z = Ang(x, y, phideg)
    Phi = atan2(y,x);
    upper = (180-phideg)*pi/180;
    lower = phideg*pi/180;
    z = double(lower<Phi & Phi<upper);
    z(Phi==upper | Phi==lower) = 0.5;
    z = z + flipud(z);

function z = Circ(x, y, D)
    r = sqrt(x.^2 + y.^2);
    z = double(r<D/2);
    z(r==D/2) = 0.5;
    
function [Eout,x,dx] = ft(Ein, z, a, u, lambda, nx,sense)

%This function finds the Fourier Transform by evaluating the integral
%instead of using FFT.  There is propogation.  This is a 2d function.
% Modified by Alexis Carlotti - Feb.17 2012

% Ein = electric field at at first plane
% z = distance of propagation in meters
% a = size of first plane (aperture) in meters
% u = size of plane 2 in meters
% lambda = wavelength
% nx = number of points in the second plane
% sense = 1 for FT, -1 for Inverse FT

nxi = length(Ein);
dxi = a/nxi;
dx = u/nx;

xi = [-(nxi-1)/2:(nxi-1)/2]*dxi;
x = [-(nx-1)/2:(nx-1)/2]*dx;
x_xi = x'*xi;

if sense~=1
    if sense ~=-1
        disp('Last parameter has to be 1 (FT), or -1 (Inverse FT); default choice is 1')
        sense=1;
    end
end

Eout=exp(2*1i*pi*z/lambda)/(1i*lambda*z)*exp(-2*sense*pi*1i*x_xi/lambda/z)*Ein*exp(-2*sense*pi*1i*x_xi'/lambda/z)*dxi*dxi;
