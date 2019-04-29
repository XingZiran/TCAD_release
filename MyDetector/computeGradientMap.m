function [LabelMap] = computeGradientMap(Ix,Iy,nBins)
    %% setup variables
    x_max = size(Ix,2);
    y_max = size(Ix,1);
    AngleVector = zeros(nBins,2);
    LabelMap = zeros(y_max,x_max);
    %% generate angles
    delta = 2 * pi / nBins;
    Angles = 0 : delta : 2 * pi - delta;
    for i = 1 : nBins
        AngleVector(i,:) = [cos(Angles(i)),sin(Angles(i))];
    end
    epson = cos(delta / 2);
    %% classify
    Gradients = [reshape(Ix,x_max * y_max,1),reshape(Iy,x_max * y_max,1)];
    % normalize each row to unit
    Gradients = Gradients ./ repmat(sqrt(sum(Gradients.^2,2)),1,size(Gradients,2));
    for iBins = 1 : nBins
        dotRes = Gradients * AngleVector(iBins,:)';
        ids = (dotRes > epson);
        LabelMap(ids) = iBins;
    end
end