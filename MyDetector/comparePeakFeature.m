function [ isExactSame,isSubSame ] = comparePeakFeature( pi,pj,nBins,dTh )
    % setup
    isExactSame = false;
    isSubSame = false;
    
    pi_halfpi = pi;
    pj_halfpi = pj;
    % since circular period, we should reform them to same region
    if pi(1) >= nBins / 4
        pi_halfpi = pi - nBins / 4 - 1;
    end
    if pj(1) >= nBins / 4
        pj_halfpi = pj - nBins / 4 - 1;
    end
    
    % get the difference of histogram
    Error1 = abs(pi - pj);
    Error2 = abs(pi_halfpi - pj);
    Error3 = abs(pi - pj_halfpi);
    
    % compare
    if norm(Error1) < dTh || norm(Error2) < dTh|| norm(Error3) < dTh
        isExactSame = true;
        return;
    elseif TripleTest(pi,pj,nBins) == true || TripleTest(pi,pj_halfpi,nBins) == true || TripleTest(pi_halfpi,pj,nBins)
        isSubSame = true;
    end
end
function [isOk] = TripleTest(Pi,Pj,nBins)
    isOk = false;
    SubTh = nBins / 18;
    Pis = [Pi([1,2,3]);Pi([1,2,4]);Pi([1,3,4]);Pi([2,3,4])];
    Pjs = [Pj([1,2,3]);Pj([1,2,4]);Pj([1,3,4]);Pj([2,3,4])];
    for i = 1 : 4
        for j = 1 : 4
            if norm(Pis(i,:) - Pjs(j,:)) < SubTh
                isOk = true;
                return;
            end
        end
    end
end
