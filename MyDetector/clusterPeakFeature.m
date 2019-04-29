function [ResIdx,MainPeakFeature] = clusterPeakFeature(PeakFeatureInEachHistogram,nBins,dTh,piTh)
    %% setup
    nPoints = size(PeakFeatureInEachHistogram,1);
    
    %% cluster peaks, i.e., cluster orientation
    Statistic = zeros(nPoints,1);
    MainMembers = cell(nPoints,1);
    SubMembers = cell(nPoints,1);
    for i = 1 : nPoints
        pi = PeakFeatureInEachHistogram(i,:);
        if min(pi) > nBins
            continue;
        end
        for j = 1 : nPoints
            pj = PeakFeatureInEachHistogram(j,:);
            if min(pj) > nBins
                continue;
            end
            [isExactSame,isSubSame] = comparePeakFeature(pi,pj,nBins,dTh);
            % filter
            if isExactSame == true
                Statistic(i) = Statistic(i) + 1;
                MainMembers{i} = [MainMembers{i},j];
            elseif isSubSame == true
                SubMembers{i} = [SubMembers{i},j];
            end
        end
    end
    %% find the biggest cluster as our result
    [~,iStatistic] = sort(Statistic,'descend');
    Imax = 0;
    for i = 1 : length(iStatistic)
        tmp = PeakFeatureInEachHistogram(iStatistic(i),:);
        piCheck1 = tmp(3) - tmp(1);
        piCheck2 = tmp(4) - tmp(2);
        if abs(piCheck1 - nBins / 2) < piTh && abs(piCheck2 - nBins / 2) < piTh
            Imax = iStatistic(i);
            break;
        end
    end
    if Imax > 0   
        % expand one step from main members
        RootIdx = MainMembers{Imax};
        tmpIdx = cell(1);
        tmpIdx{1} = RootIdx;
        for i = 1 : length(RootIdx)
            tmpIdx{1} = [tmpIdx{1},MainMembers{RootIdx(i)}];
        end
        tmpIdx{1} = [tmpIdx{1},SubMembers{Imax}];
        ResIdx = unique(tmpIdx{1});
        
        % rotate 90 degree for large perspective
        MainPeakFeature = PeakFeatureInEachHistogram(Imax,:);
        MainPeakFeature = MainPeakFeature + nBins / 4;
        for i = 1 : length(MainPeakFeature)
            if MainPeakFeature(i) > nBins
                MainPeakFeature(i) = MainPeakFeature(i) -  nBins;
            end
        end
        MainPeakFeature = sort(MainPeakFeature);
    else
        ResIdx = [];
        MainPeakFeature = [0,0,0,0];
    end
end