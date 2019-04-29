function [ Histogram ] = computeHistogramFromHAT(HAT,point,width)
    %% setup
    x_max = size(HAT,2);
    y_max = size(HAT,1);
    nBins = size(HAT,3);
    Histogram = zeros(1,nBins);
    % boundry
    x = point(1);
    y = point(2);
    min_x_coord = max(0,x - width);
    min_y_coord = max(0,y - width);
    max_x_coord = min(x_max + 1,x + width);
    max_y_coord = min(y_max + 1,y + width);
    % check validation
    if min_x_coord < 1 || min_y_coord < 1 || max_x_coord > x_max ||  max_y_coord > y_max
        return;
    end
    %% compute just like summed area table S = D + A - B - C
    tmp = HAT(max_y_coord,max_x_coord,:) + HAT(min_y_coord,min_x_coord,:) - ...
        HAT(min_y_coord,max_x_coord,:) - HAT(max_y_coord,min_x_coord,:);
    Histogram = reshape(tmp,1,nBins);
end

