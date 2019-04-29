function [ ID ] = findPtsID( I_Pts,boardSize,squareSize,Pattern,PatternPts,PatternMatrixSize,ThisImage,isDebug )
% I_Pts : image points detected from ThisImage, I_Pts is less than 
% boardSize : size of board which I_Pts lies on
% squareSize : pixel size of Pattern points grid
% Pattern : refrence pattern
% ThisImage : current input image
    %% setup
    nPts = size(I_Pts,1);
    ID = zeros(nPts,1);
    %% warp images
    [J,H] = warpImage(ThisImage,I_Pts,boardSize,squareSize);
    J_Pts = transformPointsForward(H,I_Pts);
    if isDebug
        figure;imshow(J);
        hold on;
        for i = 1 : size(J_Pts,1)
            label = sprintf('%d',i);
            text(J_Pts(i,1), J_Pts(i,2), label, 'BackgroundColor', [1 1 1]);
        end
    end
    %% get refrece area
    minY = uint32(min(J_Pts(:,2)));
    maxY = uint32(max(J_Pts(:,2)));
    minX = uint32(min(J_Pts(:,1)));
    maxX = uint32(max(J_Pts(:,1)));
%     patch = J(min(J_Pts(:,2)):max(J_Pts(:,2)),min(J_Pts(:,1)):max(J_Pts(:,1)));
    patch = J(minY:maxY,minX:maxX);
    patch = flipud(patch);% filp to demirror, since warpImage always reture the mirrored
    
     %% Use Cross-Correlation to Find Template in Image
    [patchR, patchC, ~] = size(patch); 
    [patternR, patternC, ~] = size(Pattern); 
    
    if patchR<patternR && patchC<patternC
        c = normxcorr2(patch,Pattern);
        c_max = max(c(:));

        patch180 = imrotate(patch,180);
        c180 = normxcorr2(patch180,Pattern);
        c180_max = max(c180(:));
    else
        c_max = -Inf;
        c180_max = -Inf;
    end
    
    if patchC<patternR && patchR<patternC
        patch90 = imrotate(patch,90);
        c90 = normxcorr2(patch90,Pattern);
        c90_max = max(c90(:));

        patch270 = imrotate(patch,270);
        c270 = normxcorr2(patch270,Pattern);
        c270_max = max(c270(:));
    else
        c90_max = -Inf;
        c270_max = -Inf;
    end
        
    %% find the optimal matching
    [~,max_ind]=max([c_max,c90_max,c180_max,c270_max]);
    if max_ind == 2
        patch = patch90;
        c = c90;
    elseif max_ind == 3
        patch = patch180;
        c = c180;
    elseif max_ind == 4
        patch = patch270;
        c = c270;
    end
    % Find peak in cross-correlation. 
    [ypeak, xpeak] = find(c==max(c(:)));
    %%
    yoffSet = ypeak - size(patch,1);
    xoffSet = xpeak - size(patch,2);
    %% Display matched area.
    if isDebug
        figure;
        hAx  = axes;
        imshow(Pattern,'Parent', hAx);
        imrect(hAx, [xoffSet+1, yoffSet+1, size(patch,2), size(patch,1)]);
    end
    %% find map between I_Pts and Pattern points
    % step 1 : get rectangle 4 corner coords
    Mdl = KDTreeSearcher(PatternPts);
    
    UpLeft = [xoffSet,yoffSet];
    ID_upLeft = knnsearch(Mdl,UpLeft,'K',1);
    
    UpRight = [xoffSet + size(patch,2),yoffSet];
    ID_upRight = knnsearch(Mdl,UpRight,'K',1);
    
    DownLeft = [xoffSet,yoffSet + size(patch,1)];
    ID_downLeft = knnsearch(Mdl,DownLeft,'K',1);
    
    DownRight = [xoffSet + size(patch,2),yoffSet + size(patch,1)];
    ID_downRight = knnsearch(Mdl,DownRight,'K',1);
    
    
    % step 2 : find matched grid on pattern and get its size
    [Sub_upLeft_r,Sub_upLeft_c] = ind2sub(PatternMatrixSize,ID_upLeft);
    [Sub_upRight_r,Sub_upRight_c] = ind2sub(PatternMatrixSize,ID_upRight);
    [Sub_downLeft_r,Sub_downLeft_c] = ind2sub(PatternMatrixSize,ID_downLeft);
    [Sub_downRight_r,Sub_downRight_c] = ind2sub(PatternMatrixSize,ID_downRight);
    
    if Sub_upLeft_r ~= Sub_upRight_r || Sub_downLeft_r ~= Sub_downRight_r
        return;
    else
        Up_LeftRightSize = Sub_upRight_c - Sub_upLeft_c + 1;
    end
    
    if Sub_upLeft_c ~= Sub_downLeft_c || Sub_upRight_c ~= Sub_downRight_c
        return;
    else
        Left_UpDownSize = Sub_upRight_r - Sub_downRight_r + 1;
    end
    
    % step 3 : compare grid size
    matrixSize = boardSize - 1;% matrixSize(1) for nX
    if max_ind == 1% 0 deg
        nX = matrixSize(1);
        nY = matrixSize(2);
    elseif max_ind == 2% 90 deg
        nX = matrixSize(2);
        nY = matrixSize(1);
    elseif max_ind == 3% 180 deg
        nX = matrixSize(1);
        nY = matrixSize(2);
    else % 270 deg
        nX = matrixSize(2);
        nY = matrixSize(1);
    end

    if Up_LeftRightSize ~= nY || Left_UpDownSize ~= nX
        return;%LeftRight for y-axis, UpDown for x-axis
    end
    
    % step 4 : get points ID
    PatternIndMatrix = reshape(1:PatternMatrixSize(1) * PatternMatrixSize(2),PatternMatrixSize);
    PatchIndMatrix = PatternIndMatrix(Sub_downLeft_r:Sub_upLeft_r,Sub_downLeft_c:Sub_downRight_c);
    if max_ind == 1% 0 deg
        
    elseif max_ind == 2% 90 deg
        PatchIndMatrix = fliplr(PatchIndMatrix);
        PatchIndMatrix = PatchIndMatrix';
    elseif max_ind == 3% 180 deg
        PatchIndMatrix = rot90(PatchIndMatrix,2);
    else % 270 deg
        PatchIndMatrix = flipud(PatchIndMatrix);
        PatchIndMatrix = PatchIndMatrix';
    end
    ID = PatchIndMatrix(:);
end