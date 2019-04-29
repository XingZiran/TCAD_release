function [ loc_filtered ] = removeInvalidPeaks( metric,loc,patchsize,width,ratio)
% metric : data matrix m_x_n
% threshold : threshold in peak searching
% loc : result of peak searching result, in pixel coordinate
% width : maximum width threshold of decreasing ratio
% ratio : decreasing ratio larger than 1
% patchsize : size of maximum patch, define a (2 x patchsize + 1)^2 window
    %% setup variables
    nPts = size(loc,1);
    loc_filtered = zeros(nPts,2);
    
    %% test each point to know if this point is valid
    k = 0;% size of valid peaks
    for i = 1 : nPts
        % test local maximum
        IsPeakValid = PeakMaximumTest(metric,loc,i,patchsize);
        IsWidthValid = PeakWidthTest(metric,loc,i,width,ratio);
        % remove invalid Peaks
        if IsWidthValid == true && IsPeakValid == true 
            k = k + 1;
            loc_filtered(k,:) = loc(i,:); 
        end
    end
    %% generate output
    if k > 0
        loc_filtered = loc_filtered(1:k,:);
    else
        loc_filtered = [];
    end
end
%% 
function [IsValid] = PeakMaximumTest(metric,loc,i,width)
    [m,n] = size(metric);% m for #rows and n for #col
    % change the boundry to image cordinate
    x_max = n;
    y_max = m;
    % get 8-neibors defined by width of window
    x = loc(i,1);
    y = loc(i,2);
    % find up-bound
    up = max(y - width,1);
    down = min(y + width,y_max);
    left = max(x - width,1);
    right = min(x + width,x_max);
    % get region maximum
    maximum = max(max(metric(up:down,left:right)));
    Peak = metric(y,x);
    if Peak >= maximum
        IsValid = true;
    else
        IsValid = false;
    end
end
%%
function [IsValid] = PeakWidthTest(metric,loc,i,width,ratio)
    [m,n] = size(metric);% m for #rows and n for #col
    % change the boundry to image cordinate
    x_max = n;
    y_max = m;
    % get 16-neibors defined by width of window
    x = loc(i,1);
    y = loc(i,2);
    neiborLoc = generateCoordinate(x_max,y_max,x,y,width);
    % test if every 8-neibors has less than half intensity of peak
    Peak = metric(y,x);
    IsValid = true;
    for j = 1 : 8
        xj = neiborLoc(j,1);
        yj = neiborLoc(j,2);
        % only consider those neibor points inside the region
        if xj > 0 && yj > 0
            neibor = metric(yj,xj);
            if neibor > Peak / ratio
                % this point is not valid
                IsValid = false;
                break;
            end
        end
    end
end
%%
function [neiborLoc] = generateCoordinate(x_max,y_max,x,y,width)
    neiborLoc = zeros(8,2);
    % up-left
    tmp_X = x - width;
    tmp_Y = y - width;
    if tmp_X > 0 && tmp_Y > 0
        neiborLoc(1,1) = tmp_X;
        neiborLoc(1,2) = tmp_Y;
    end
    % left
    tmp_X = x - width;
    tmp_Y = y;
    if tmp_X > 0
        neiborLoc(2,1) = tmp_X;
        neiborLoc(2,2) = tmp_Y;
    end
    % down-left
    tmp_X = x - width;
    tmp_Y = y + width;
    if tmp_X > 0 && tmp_Y <= y_max
        neiborLoc(3,1) = tmp_X;
        neiborLoc(3,2) = tmp_Y;
    end
    % up
    tmp_X = x;
    tmp_Y = y - width;
    if tmp_Y > 0
        neiborLoc(4,1) = tmp_X;
        neiborLoc(4,2) = tmp_Y;
    end
    % down
    tmp_X = x;
    tmp_Y = y + width;
    if tmp_Y <= y_max
        neiborLoc(5,1) = tmp_X;
        neiborLoc(5,2) = tmp_Y;
    end
    % up-right
    tmp_X = x + width;
    tmp_Y = y - width;
    if tmp_X <= x_max && tmp_Y > 0
        neiborLoc(6,1) = tmp_X;
        neiborLoc(6,2) = tmp_Y;
    end
    % right
    tmp_X = x + width;
    tmp_Y = y;
    if tmp_X <= x_max
        neiborLoc(7,1) = tmp_X;
        neiborLoc(7,2) = tmp_Y;
    end
    % down-right
    tmp_X = x + width;
    tmp_Y = y + width;
    if tmp_X <= x_max && tmp_Y <= y_max
        neiborLoc(8,1) = tmp_X;
        neiborLoc(8,2) = tmp_Y;
    end
end