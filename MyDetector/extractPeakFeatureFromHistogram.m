function [PeakFeature] = extractPeakFeatureFromHistogram(Histogram)
    % setup
    nBins = length(Histogram);
    PeakFeature = [1000 1000 1000 1000];
    OutstandingRiatio = 0.1;
    % make CircularHistogram circular since the Histogram means sample on unit circle
    CircularHistogram = [Histogram(nBins),Histogram,Histogram(1)];
    [pks.value,pks.locs,pks.width,pks.prominence] = findpeaks(CircularHistogram);
    pks.locs = pks.locs - 1;
    %% select 4 maximum peaks
    [Res,Idx] = sort(pks.value,'descend');
    if length(Idx) < 4
        return;
    elseif length(Idx) >= 5
        if (Res(3) - Res(5)) < OutstandingRiatio * (Res(1) - Res(5))
            return;
        end
    end
    PeakFeature = sort(pks.locs(Idx(1:4)));
end