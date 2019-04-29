function [PeakFeature] = filterPtsByHistogram( HAT,pointCoords,nRefWidth,step,dth)
    %% compute histogram at different size
    Histograms_step1 = computeHistogramFromHAT(HAT,pointCoords,nRefWidth - step);
    Histograms_step2 = computeHistogramFromHAT(HAT,pointCoords,nRefWidth);
    Histograms_step3 = computeHistogramFromHAT(HAT,pointCoords,nRefWidth + step);
    nBins = length(Histograms_step2);
    %% Extract peak features
    PeakFeature_minus = extractPeakFeatureFromHistogram(Histograms_step1);
    PeakFeature = extractPeakFeatureFromHistogram(Histograms_step2);
    PeakFeature_plus = extractPeakFeatureFromHistogram(Histograms_step3);
    %% compare
    [isStep1Same,isStep1SubSame] = comparePeakFeature(PeakFeature,PeakFeature_minus,nBins,dth);
    [isStep3Same,isStep3SubSame] = comparePeakFeature(PeakFeature,PeakFeature_plus,nBins,dth);
    if ((isStep1Same || isStep1SubSame) || (isStep3Same || isStep3SubSame)) == false
        PeakFeature = [1000,1000,1000,1000];
    end
end