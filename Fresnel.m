function [Eout,x,dx] = Fresnel(Ein, z, a, u, lambda, nx, trigger)

% Fresnel propagation ; You can choose a classic Fresnel propagation or 
% a smarter move (FT of Ein, propagation in the Fourier plane, then FT-1) 
% Modified by Alexis Carlotti - Feb.17 2012
% Modified by Aaron Lemmer - May 17 2012

% Ein = electric field at at first plane
% z = distance of propagation in meters
% a = size of plane 1 in meters
% u = size of plane 2 in meters
% lambda = wavelength
% nx = number of points in the second plane
% trigger = 0 for simple Fresnel transform
% trigger !=0 for angular spectrum

if trigger==0

    nxi = length(Ein);  %number of elements in input matrix, Ein (electric
                        %field at input plane), taken along longest
                        %dimension of Ein

        
    dxi = a/nxi;   %step size in input plane (plane 1)
    dx = u/nx;   %step size in output plane (plane 2)
    xi = [-nxi/2+0.5:nxi/2-0.5]*dxi;   %creates x vector in input plane,
                                       %using step size defined above
    x = [-nx/2+0.5:nx/2-0.5]*dx;   %creates x vector in output plane
    x_xi = x'*xi;   %

    
    [XI,ETA] = meshgrid(xi);   %creates square (xi by xi) grid of (xi,yi)
                             %coordinate pairs in plane 1
    [X,Y] = meshgrid(x);   %creates square (x by x) grid of (x,y)
                           %coordinate pairs in plane 2

    PropFactor = exp(1i*pi*(XI.^2 + ETA.^2)/(lambda*z));
    
    %figure(100)
    %imagesc(xi, xi, cos(angle(Ein) + pi*(XI.^2 + ETA.^2)/(lambda*z)))
    %colormap('gray');
    %axis square
    %colorbar

    Eout = exp(-2*pi*1i*x_xi/lambda/z)*(Ein.*PropFactor)*exp(-2*pi*1i*x_xi'/lambda/z)*dxi*dxi;

    Eout = exp(2*1i*pi*z/lambda)/(1i*lambda*z)*exp(1i*pi*(X.^2 + Y.^2)/(lambda/z)).*Eout;

else
    
    Nring = 100;     % Number of airy rings in the Fourier plane
    Npix = 10;       % Number of pixels per ring
    
    Dint = 2*Nring*lambda*z/a; % Width of the Fourier plane
    Nint = 2*Npix*Nring;       % Number of pixels in the Fourier plane

    [four,x_int,dx_int] = ft(Ein,z,a,Dint,lambda,Nint,1);
    
    [X_int,Y_int] = meshgrid(x_int);
    
    angular = exp(-1i*pi*z*lambda*(X_int.^2 + Y_int.^2));
    
    %figure(100)
    %imagesc(x_int, x_int, cos(angle(four) + pi*z*lambda*(X_int.^2 + Y_int.^2)))
    %axis square
    %colorbar

    [Eout,x,dx] = ft(angular.*four,z,Dint,u,lambda,nx,-1);

    Eout = exp(2*1i*pi*z/lambda)*Eout;
    
end
