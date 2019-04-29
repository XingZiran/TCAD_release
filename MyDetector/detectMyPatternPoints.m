function [ points, boardSize ] = detectMyPatternPoints( Igray,isDebug )
%DE Summary of this function goes here
%   Detailed explanation goes here
    I = im2single(Igray);
    [m,n] = size(I);

    % Bandwidth of the gaussian filter for corner detection
    % If a checkerboard is not detected in a high-resolution image, increase
    % the value of sigma
    sigma = 2; 

    minCornerMetric = 0.15; % threshold for corner metric

    
    ParaSet.nBins = 64;             % number of bins for orientation estimation
    ParaSet.width = 2;
    ParaSet.ratio = 1.5;            % for width test
    if m < 800 || n < 800
        ParaSet.patchsize = 3;         % for large region maximum test
        ParaSet.nLongPath = 10;
        ParaSet.nRefWidth = 10;         % for histogram computation
    else
        ParaSet.patchsize = 10;         % for large region maximum test
        ParaSet.nLongPath = 10;
        ParaSet.nRefWidth = 10;         % for histogram computation
    end
    
    [points, boardSize] = detectCheckerboard(I, sigma, minCornerMetric,ParaSet,isDebug);

end

