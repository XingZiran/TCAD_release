function [ HAT ] = computeHistogramAreaTableFast( LabelMap,nBins )
    m = size(LabelMap,1);
    n = size(LabelMap,2);
    HAT = zeros(m,n,nBins);
    %% initial : build the histogram map
    indices = find(LabelMap > 0);
    [row,col] = ind2sub([m,n],indices);
    for i = 1 : length(indices)
        HAT(row(i),col(i),LabelMap(indices(i))) = 1;
    end
    HAT = cumsum(cumsum(HAT,1),2);
end