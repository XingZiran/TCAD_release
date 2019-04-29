function [cxy, c45, Ix, Iy] = secondDerivCornerMetric(I, sigma)
%#codegen
assert(ismatrix(I));

% Low-pass filter the image
G = fspecial('gaussian', round(sigma * 7)+1, sigma);
Ig = imfilter(I, G, 'conv');

derivFilter = [-1 0 1];

% first derivatives
Iy = imfilter(Ig, derivFilter', 'conv');
Ix = imfilter(Ig, derivFilter, 'conv');
    
% first derivative at 45 degrees
I_45 = Ix * cos(pi/4) + Iy * sin(pi/4);
I_n45 = Ix * cos(-pi/4) + Iy * sin(-pi/4);

% second derivative
Ixy = imfilter(Ix, derivFilter', 'conv');

I_45_x = imfilter(I_45, derivFilter, 'conv');
I_45_y = imfilter(I_45, derivFilter', 'conv');    

I_45_45 = I_45_x * cos(-pi/4) + I_45_y * sin(-pi/4);

% suppress the outer corners
cxy = sigma^2 * abs(Ixy) - sigma * (abs(I_45) + abs(I_n45));
cxy(cxy < 0) = 0;
c45 = sigma^2 * abs(I_45_45) - sigma * (abs(Ix) + abs(Iy));
c45(c45 < 0) = 0;


