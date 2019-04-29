function [ filtedCorners ] = filterPtsByLongPath(CornerMetric,CornerCandidates,refWidth,isDebug)
% filter corner candidates by long path on metric data
    %% make binary Metric
    biMetic = imbinarize(CornerMetric);
    if isDebug
        figure;imshow(biMetic);
        title('Binarize reult!');
    end
    cc = bwconncomp(biMetic,8);% 8-connectivity
    %% remove large cc by its bounding box
    stats = regionprops(cc,'BoundingBox');
    nRegions = length(cc.PixelIdxList);
    filtedRegions = zeros(nRegions,4);
    k = 0;
    for i = 1 : nRegions
        BBX = stats(i).BoundingBox;
        x_width = BBX(3);
        y_width = BBX(4);
        if x_width <= refWidth && y_width <= refWidth %% check if size meet require
            k = k + 1;
            filtedRegions(k,:) = BBX;
        end
    end
    %% filte candidates by check if each point is inside the region and belongs to its cc
    filtedCorners = CornerCandidates;
    nCorners = 0;
    for i = 1 : size(CornerCandidates,1)
        x = CornerCandidates(i,1);% image pixel is different ???
        y = CornerCandidates(i,2);
        for j = 1 : k
            BBX = filtedRegions(j,:);
            ulx = BBX(1);
            uly = BBX(2);
            x_width = BBX(3);
            y_width = BBX(4);
            if x >= ulx && x <= ulx + x_width && y >= uly && y <= uly + y_width
                nCorners = nCorners + 1;
                filtedCorners(nCorners,:) = CornerCandidates(i,:);
            end
        end
    end
    filtedCorners = filtedCorners(1:nCorners,:);
end