function [points, boardSize] = detectCheckerboard(I, sigma, peakThreshold,ParameterSet,isDebug)
    %% common process
    I = preprocessImage(I);
    [cxy, c45, Ix, Iy] = secondDerivCornerMetric(I, sigma);
    nBins = ParameterSet.nBins;         % number of bins for orientation estimation
    width = ParameterSet.width;         % for width test
    ratio = ParameterSet.ratio;         % for width test
    patchsize = ParameterSet.patchsize; % for large region maximum test
    nLongPath = ParameterSet.nLongPath; % 
    nRefWidth = ParameterSet.nRefWidth; % for histogram computation
    nStep = 2;                          % for resize histogram filtering
    dTh = nBins / 16;                   % for cluter points by histogram
    piTh = 10 / 360 * nBins;            % for remove points which their orientation is not straight
    %% detect corners in cxy and c45
    points0 = findPeaksMatlab(cxy, peakThreshold);
    
    points45 = findPeaksMatlab(c45, peakThreshold);
    
    isSkip0 = false;
    isSkip45 = false;
    
    if isempty(points0)
        isSkip0 = true;
    elseif isDebug
        figure;imshow(cxy);hold on;plot(points0(:,1),points0(:,2),'ro');
        title('Original peak detection by MATLAB, 0 deg');
    end 
    
    if isempty(points45)
        isSkip45 = true;
    elseif isDebug
        figure;imshow(c45);hold on;plot(points45(:,1),points45(:,2),'ro');
        title('Original peak detection by MATLAB, 45 deg');
    end
    %% step 1 : preprocess, remove by width and peak
    % debug
    points0 = filterPtsByLongPath(cxy,points0,nLongPath,isDebug);
    if isDebug && ~isempty(points0)
        figure;imshow(cxy);hold on;plot(points0(:,1),points0(:,2),'ro');
        title('Filtering by longest path, 0 deg');
    end 
    points0_original = points0;
    points0 = removeInvalidPeaks(cxy,points0,patchsize,width,ratio);
    if isempty(points0)
        isSkip0 = true;
    elseif isDebug
        figure;imshow(cxy);hold on;plot(points0(:,1),points0(:,2),'ro');
        title('My first suppress, 0 deg');
    end 
    
    points45 = filterPtsByLongPath(c45,points45,nLongPath,isDebug);
    if isDebug && ~isempty(points45)
        figure;imshow(cxy);hold on;plot(points45(:,1),points45(:,2),'ro');
        title('Filtering by longest path, 45 deg');
    end
    points45_original = points45;
    points45 = removeInvalidPeaks(c45,points45,patchsize,width,ratio);
    
    if isempty(points45)
        isSkip45 = true;
    elseif isDebug
        figure;imshow(c45);hold on;plot(points45(:,1),points45(:,2),'ro');
        title('My first suppress,45 deg');
    end
    %% preprocess for histogram
    if isSkip0 == true && isSkip45 == true
        points = [];
        boardSize = [];
        return;
    else
        GradientMap = computeGradientMap(Ix,Iy,nBins);
        HAT = computeHistogramAreaTableFast(GradientMap,nBins);
%         HAT = computeHistogramAreaTable(GradientMap,nBins);
    end 
    %% step 2 : process filter result which angle is 0 and 45
    % 0
    if isSkip0 == false
        PeakFeatures0 = zeros(size(points0,1),4);
        for i = 1 : size(points0,1)
            PeakFeatures0(i,:) = filterPtsByHistogram(HAT,points0(i,:),nRefWidth,nStep,dTh);
        end

        % filte the candidates peaks by clustering histogram
        [filtedIdx0,MainPeakFeature0] = clusterPeakFeature(PeakFeatures0,nBins,dTh,piTh);
        points0 = points0(filtedIdx0,:);
        
        if isDebug
            figure;imshow(cxy);hold on;plot(points0(:,1),points0(:,2),'ro');
            title('Histogram filtering result,0 deg');
        end

        % recover the checker board structure
        scores0 = cxy(sub2ind(size(cxy), points0(:, 2), points0(:, 1)));
        board0 = growCheckerboardX(points0,scores0,MainPeakFeature0,nBins,cxy);
    else
        board0 = struct('BoardIdx', zeros(3), 'BoardCoords', zeros(3,3,3), ...
                'Energy', Inf, 'isValid', false);
    end
    % 45
    if isSkip45 == false
        PeakFeatures45 = zeros(size(points45,1),4);
        for i = 1 : size(points45,1)
            PeakFeatures45(i,:) = filterPtsByHistogram(HAT,points45(i,:),nRefWidth,nStep,dTh);
        end

        % filte the candidates peaks in clustering histogram
        [filtedIdx45,MainPeakFeature45] = clusterPeakFeature(PeakFeatures45,nBins,dTh,piTh);
        points45 = points45(filtedIdx45,:);
        
        if isDebug
            figure;imshow(c45);hold on;plot(points45(:,1),points45(:,2),'ro');
            title('Histogram filtering result,45 deg');
        end
        
        % grow checkerboard by structure information
        scores45 = c45(sub2ind(size(c45), points45(:, 2), points45(:, 1)));
        board45 = growCheckerboardX(points45,scores45,MainPeakFeature45,nBins,c45);
    else
        board45 = struct('BoardIdx', zeros(3), 'BoardCoords', zeros(3,3,3), ...
                'Energy', Inf, 'isValid', false);
    end
    %% check the validation of board
    [isValid0,points0,boardSize0,boardEnergy0] = gridCheck(board0,cxy,points0_original);%
    [isValid45,points45,boardSize45,boardEnergy45] = gridCheck(board45,c45,points45_original);% 
    %% selet the best board
    points = [];
    boardSize = [0 0];
    if isValid0 && boardEnergy0 < boardEnergy45
        points = points0;
        boardSize = boardSize0;
        if isDebug
            disp('0 deg');
            figure(102);imshow(I);hold on;plot(points0(:,1),points0(:,2),'ro');hold on;
            plot(points(:,1),points(:,2),'b*');
        end
    elseif isValid45
        points = points45;
        boardSize = boardSize45;
        if isDebug
            disp('45 deg');
            figure(102);imshow(I);hold on;plot(points45(:,1),points45(:,2),'ro');hold on;
            plot(points(:,1),points(:,2),'b*');
        end
    end
end
%% TPAMI'07 estimate fundamental matrix with radial distortion and refine the points location
function [isValid,points,boardSize,boardEnergy] = gridCheck(board,metric,allCandidates)
    % setup
    isValid = false;
    points = [];
    boardSize = [0,0];
    boardEnergy = inf;
    if board.isValid == false
        return;
    end
    isValid = true;
    % get points
    [points, boardSize] = toPoints(board);
    points = subPixelLocation(metric, points);
    % generate word points
    squareSize = 500;
    worldPoints = generateCheckerboardPoints(boardSize,squareSize);
    % fit fundamental matrix and transform image points to ideal grid
    MovingPts = points;
    FixedPts = worldPoints;
    inlierTh = 0.1;
    
    % esitmate fundamental matirx, for more robust, try 2 times
    p1 = [FixedPts,ones(size(FixedPts,1),1)];
    p2 = [MovingPts,ones(size(MovingPts,1),1)];
    [fMatrix,~,Status] = estimateFundamentalMatrix(FixedPts,MovingPts,...
        'Method','LTS',...
        'InlierPercentage',91,...
        'DistanceType','Algebraic',...
        'ReportRuntimeError',false);
    if Status == 0
        d = abs(diag(p2 * fMatrix * p1'));
        d_mean = mean(d);
    else
        d_mean = inf;
    end
    
    % check status and compute result
    
    if d_mean < inlierTh
        % refine coords location if we can
        kdTree = KDTreeSearcher(allCandidates);
        for i = 1 : length(d)
            if d(i) >= 10 * d_mean
                disp('Outliner exist Start Refineing...');
                ptCoord = refineLocation(kdTree,points,i,boardSize,fMatrix,d,worldPoints(i,:));
                points(i,:) = subPixelLocation(metric, ptCoord);
            end
        end
        % 
        boardEnergy = getBoardEnergy(points,boardSize);
        if size(points,1) < 12 && boardEnergy > -7.9
            isValid = false;
            disp('Poor board...');
        end
    else
        disp('Too many outliers exist...');
        disp(d_mean);
        isValid = false;
    end
end
%% get board energy
function [boardEnergy] = getBoardEnergy(points,boardSize)
    % setup
    gridSize = boardSize - 1;
    boardEnergy = inf;
    m = gridSize(1);
    n = gridSize(2);
    if size(points,1) ~= m * n
        return;
    end
    % compute maxtriple
    maxTriple = 0;
    for i = 2 : m - 1
        for j = 2 : n - 1
            thisID = sub2ind(gridSize,i,j);
            upID = sub2ind(gridSize,i - 1,j);
            downID = sub2ind(gridSize,i + 1,j);
            leftID = sub2ind(gridSize,i,j - 1);
            rightID = sub2ind(gridSize,i,j + 1);
            Triple1 = norm(points(upID,:) + points(downID,:) - 2 * points(thisID,:))...
                ./ norm(points(upID,:) - points(downID,:));
            Triple2 = norm(points(leftID,:) + points(rightID,:) - 2 * points(thisID,:))...
                ./ norm(points(leftID,:) - points(rightID,:));
            maxTriple = max([maxTriple,Triple1,Triple2]);
        end
    end
    %% 
    boardEnergy = -1 * m * n + maxTriple * m * n;
end
%% refine the points location
function [ptCoord] = refineLocation(kdTree,pts,i,boardSize,fMatrix,d,wpts)
    % setup
    ptCoord = pts(i,:);
    d_mean = mean(d);
    % get center and radius
    neiborIDs = getNeiborID(boardSize,i);
    validNeribor = zeros(length(neiborIDs),1);
    validCoords = zeros(length(neiborIDs),2);
    for i = 1 : length(neiborIDs)
        if neiborIDs(i) > 0
            if d(neiborIDs(i)) < 10 * d_mean
                validNeribor(i) = 1;
                validCoords(i,:) = pts(neiborIDs(i),:);
            end
        end
    end
    if sum(validNeribor) < 2
        return;
    elseif sum(validNeribor) > 3
        Center = sum(validCoords,1)./sum(validNeribor);
        radius = inf;
        for i = 1 : length(neiborIDs)
            if validNeribor(i) > 0
                radius = min(radius,norm(Center - validCoords(i,:)));
            end
        end
    elseif sum(validNeribor) == 3
        if validNeribor(1) == 1 && validNeribor(3) == 1
            % up and down
            Center = (validCoords(1,:) + validCoords(3,:)) / 2;
            radius = min(norm(Center - validCoords(1,:)),norm(Center - validCoords(3,:)));
        elseif validNeribor(2) == 1 && validNeribor(4) == 1
            % left and right
            Center = (validCoords(2,:) + validCoords(4,:)) / 2;
            radius = min(norm(Center - validCoords(2,:)),norm(Center - validCoords(4,:)));
        else
            return;
        end
    else % here is 2 valid neibor condition
        if validNeribor(1) == 1 && validNeribor(3) == 1
            % up and down
            Center = sum(validCoords,1)./sum(validNeribor);
            radius = min(norm(Center - validCoords(1,:)),norm(Center - validCoords(3,:)));
        elseif validNeribor(2) == 1 && validNeribor(4) == 1
            % left and right
            Center = sum(validCoords,1)./sum(validNeribor);
            radius = min(norm(Center - validCoords(2,:)),norm(Center - validCoords(4,:)));
        elseif validNeribor(1) == 1 && validNeribor(2) == 1
            % up and left
            idUpLeft = neiborIDs(2) - 1;
            if d(idUpLeft) < 10 * d_mean
                ptUpLeft = pts(idUpLeft,:);
                Center = validCoords(1,:) + validCoords(2,:) - 2 * ptUpLeft + ptUpLeft;
                radius = min(norm(Center - validCoords(1,:)),norm(Center - validCoords(2,:)));
            else
                return;
            end
        elseif validNeribor(2) == 1 && validNeribor(3) == 1
            % left and down
            idLeftDown = neiborIDs(2) + 1;
            if d(idLeftDown) < 10 * d_mean
                ptLeftDown = pts(idLeftDown,:);
                Center = validCoords(2,:) + validCoords(3,:) - 2 * ptLeftDown + ptLeftDown;
                radius = min(norm(Center - validCoords(2,:)),norm(Center - validCoords(3,:)));
            else
                return;
            end
        elseif validNeribor(3) == 1 && validNeribor(4) == 1
            % right and down
            idRightDown = neiborIDs(4) + 1;
            if d(idRightDown) < 10 * d_mean
                ptRightDown = pts(idRightDown,:);
                Center = validCoords(3,:) + validCoords(4,:) - 2 * ptRightDown + ptRightDown;
                radius = min(norm(Center - validCoords(3,:)),norm(Center - validCoords(4,:)));
            else
                return;
            end
        elseif validNeribor(1) == 1 && validNeribor(4) == 1
            % right and up
            idRightUp = neiborIDs(4) - 1;
            if d(idRightUp) < 10 * d_mean
                ptRightUp = pts(idRightUp,:);
                Center = validCoords(1,:) + validCoords(4,:) - 2 * ptRightUp + ptRightUp;
                radius = min(norm(Center - validCoords(1,:)),norm(Center - validCoords(4,:)));
            else
                return;
            end
        else
            return;
        end
    end
    % compute new location
    IdxCell = rangesearch(kdTree,Center,radius * 0.9);
    Indices = IdxCell{1};
    opt = inf;
    optIdx = 0;
    x = kdTree.X;
    for i = 1 : length(Indices)
        if opt > abs([x(Indices(i),:),1] * fMatrix * [wpts,1]') + 0.005
            optIdx = Indices(i);
            opt = abs([x(Indices(i),:),1] * fMatrix * [wpts,1]');
        end
    end
    ptCoord = x(optIdx,:);
end
%% find neibor ID in board, 4-neibor
function [neiborIDs] = getNeiborID(boardSize,i)
    neiborIDs = zeros(4,1);
    matrixSize = boardSize - 1;
    [r,c] = ind2sub(matrixSize,i);
    if r > 1
        neiborIDs(1) = sub2ind(matrixSize,r - 1,c);%up
    end
    if c > 1
        neiborIDs(2) = sub2ind(matrixSize,r,c - 1);%left
    end
    if r < matrixSize(1)
        neiborIDs(3) = sub2ind(matrixSize,r + 1,c);%down
    end
    if c < matrixSize(2)
        neiborIDs(4) = sub2ind(matrixSize,r,c + 1);%right
    end
end
%% --------------------------------------------------------------------------
function I = preprocessImage(I)

% add a tiny amount of noise in case the image is "perfect"
% if the squares are prefectly flat, then local maxima detection may not
% work correctly.
I = I + rand(size(I)) * 1e-10;

end
%% --------------------------------------------------------------------------
function [points, boardSize] = toPoints(this)
% returns the points as an Mx2 matrix of x,y coordinates, and
% the size of the board

    if any(this.BoardIdx(:) == 0)
        points = [];
        boardSize = [0 0];
        return;
    end

    numPoints = size(this.BoardCoords, 1) * size(this.BoardCoords, 2);
    points = zeros(numPoints, 2);
    x = this.BoardCoords(:, :, 1)';
    points(:, 1) = x(:);
    y = this.BoardCoords(:, :, 2)';
    points(:, 2) = y(:);
    boardSize = [size(this.BoardCoords, 2)+1, size(this.BoardCoords, 1)+1];
end