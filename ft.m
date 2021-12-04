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