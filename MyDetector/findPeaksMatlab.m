function loc = findPeaksMatlab(metric,quality)
% Returns local maxima in the image metric. quality is a threshold
% expressed as a fraction of the maximum value of metric.

%#codegen
maxMetric = max(metric(:));
if maxMetric <= eps(0)
    loc = zeros(0,2, 'single');
else
    
    bw = imregionalmax(metric, 8);
    
    threshold = quality * maxMetric;
    bw(metric < threshold) = 0;
    bw = bwmorph(bw, 'shrink', Inf);
    
    % Exclude points on the border
    bw(1, :) = 0;
    bw(end, :) = 0;
    bw(:, 1) = 0;
    bw(:, end) = 0;
    
    % Find location of the peaks
    idx = find(bw);
    loc = zeros([length(idx) 2], 'like', metric );
    [loc(:, 2), loc(:, 1)] = ind2sub(size(metric), idx);
end