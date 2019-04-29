function [ IPts,ID,boardSize ] = detectAndLabelBoard(I,Pattern,PatternPts,PatternMatrixSize,squareSize,...
    IsNeedID,IsDebug,IsVisualization)
    %% setup
    if IsVisualization
        figure;imshow(Pattern);hold on;
        for i = 1 : size(PatternPts,1)
            label = sprintf('%d',i);
            text(PatternPts(i,1), PatternPts(i,2), label,'BackgroundColor', [1 1 1]);
        end 
    end
    %% detect and label saddle points on I
    [IPts,boardSize] = detectMyPatternPoints(I,IsDebug);
    if ~isempty(IPts)
        if IsNeedID
            ID = findPtsID(IPts,boardSize,squareSize,Pattern,PatternPts,PatternMatrixSize,I,IsDebug);
        else
            ID = zeros(size(IPts,1),1);
        end
        if IsVisualization
            figure;imshow(Image);hold on;
            for i = 1 : size(IPts,1)
                label = sprintf('%d',ID(i));
                text(IPts(i,1), IPts(i,2), label,'BackgroundColor', [1 1 1]);
            end
            disp('Detection complete!');
        end
    else
        ID = [];
        if IsVisualization
            disp('No checkerboard detected!');
        end
    end
end