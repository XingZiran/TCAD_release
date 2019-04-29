function [J,H] = warpImage(Image,I_Pts,boardSize,squareSize)
% Input : Image , m_x_n image to be warp
%       : I_Pts , Checkerboard points which to be warp nPts_x_2
%       : boardSize, boardSize for detected checkerboard
%       : squareSize, nPixels of model pattern square
    %% generate checkerboard
    d = 100;
    [worldPoints] = generateCheckerboardPoints(boardSize,squareSize);
    worldPoints = worldPoints + d;
    %% compute Homograph
    MovingPts = I_Pts;
    FixedPts = worldPoints;
    H = fitgeotrans(MovingPts, FixedPts, 'projective');% moving relative to fix
    %% add 
%     addpath('./RANSAC_Homography/');
%     [Hr, inliers] = ransacfithomography(MovingPts', FixedPts', 0.1);
%     Hr = Hr' ./ Hr(3,3);
%     H.T = Hr;
    roix = worldPoints(size(worldPoints,1),1) + d;
    roiy = worldPoints(size(worldPoints,1),2) + d;
    R = imref2d([roiy,roix]);
    J = imwarp(Image,H,'OutputView',R);
end