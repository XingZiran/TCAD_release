function [ Histogram ] = computeGradientHistogram( Ix,Iy,x,y,Width,NumberOfBins )
% Ix : dI/dx matrix with same size if original images
% Iy : dI/dy matrix with same size if original images
% Width : the half width reference window
% NumberOfBins : number of samples in angle
    %% setup variables
    x_max = size(Ix,2);
    y_max = size(Ix,1);
    Histogram = zeros(1,NumberOfBins);
    AngleVector = zeros(NumberOfBins,2);
    min_x_coord = max(0,x - Width);
    min_y_coord = max(0,y - Width);
    max_x_coord = min(x_max + 1,x + Width);
    max_y_coord = min(y_max + 1,y + Width);
    % check validation
    if min_x_coord < 1 || min_y_coord < 1 || max_x_coord > x_max ||  max_y_coord > y_max
        return;
    end
    %% generate angles
    delta = 2 * pi / NumberOfBins;
    Angles = 0 : delta : 2 * pi - delta;
    for i = 1 : NumberOfBins
        AngleVector(i,:) = [cos(Angles(i)),sin(Angles(i))];
    end
    epson = cos(delta / 2);
    %% Test each gradient of each pixel
    for x = min_x_coord : 1 : max_x_coord
        for y = min_y_coord : 1 : max_y_coord
            Current_gradient = [Ix(y,x),Iy(y,x)] ./ norm([Ix(y,x),Iy(y,x)]);
            for i = 1 : NumberOfBins
                if dot(Current_gradient, AngleVector(i,:)) > epson
                    Histogram(i) = Histogram(i) + 1;
                    break;
                end
            end
        end
    end
end