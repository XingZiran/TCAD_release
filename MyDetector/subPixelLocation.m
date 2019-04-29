function loc = subPixelLocation(metric, loc)

%#codegen 

for id = 1: size(loc,1)
    loc(id,:) = subPixelLocationImpl(metric, loc(id,:));
end

function subPixelLoc = subPixelLocationImpl(metric, loc)
% subPixelLocation(metric, loc) fits a quadratic function to a patch in the
% image metric around the pixel specified by loc. metric is a grayscale
% 2-D image, and loc is a two-element vector [x y]. subPixelLoc is a 
% two-element vector [x y], corresponding to the maximum of the quadratic 
% function.

% The size of the patch is halfPatchSize * 2 + 1
halfPatchSize = 2;

% Check if the patch is outside the image
if any(loc(:) < halfPatchSize + 1) || loc(1) > size(metric, 2) - halfPatchSize - 1 ...
    || loc(2) > size(metric, 1) - halfPatchSize -1
    subPixelLoc = single(loc);
    return;
end

% Get the patch
patch = metric(loc(2)-halfPatchSize:loc(2)+halfPatchSize, ...
               loc(1)-halfPatchSize:loc(1)+halfPatchSize);

% Create the matrix for the normal equations used to solve for the 
% coefficients of the quadratic function. This matrix depends only on the 
% patch size, and thus needs to be computed only once.
persistent X;
if isempty(X)
    X = createNormalEquationsMatrix(halfPatchSize);
end

% Get a least-squares solution for the coefficients
y = patch(:);
beta = X * y;

% f(x,y) = Ax^2 + By^2 + Cx + Dy + Exy + F
A = beta(1);
B = beta(2);
C = beta(3);
D = beta(4);
E = beta(5);
% F = beta(6), but we do not need it

% Solve for the maximum of the quadratic, in patch-based coordinates
x = -(2*B*C - D*E) / (4*A*B - E^2);
y = -(2*A*D - C*E) / (4*A*B - E^2);
if ~isfinite(x) || abs(x) > halfPatchSize || ~isfinite(y) || abs(y) > halfPatchSize
    x = single(0);
    y = single(0);
end

% Get the sub-pixel location
subPixelLoc = single(loc) + [x y];
 
function X = createNormalEquationsMatrix(halfPatchSize)
% X = [1 1 -1 -1  1 1;
%      0 1  0 -1  0 1;
%      1 1  1 -1 -1 1; % end row 1
%      1 0 -1  0  0 1;
%      1 0  1  0  0 1; 
%      0 0  0  0  0 1; % end row 2
%      1 1 -1  1 -1 1;
%      0 1  0  1  0 1;
%      1 1  1  1  1 1];

v = -halfPatchSize : 1 : halfPatchSize;
[x, y] = meshgrid(v, v);
x = x(:);
y = y(:);
X = [x.^2, y.^2, x, y, x.*y, ones(size(x))];
X = (X' * X) \ X';
     
