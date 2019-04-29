function board = growCheckerboardX(points,scores,MainPeakFeature,nBins,I)
% points : input corner candidates nPts_x_2
% score : Response value of each corner candidates nPts_x_1
% MainPeakFeature : nPts_x_4 peak feature of each corner candidates'histogram
% nBins : number of bins in each histogram
% I : only need for debug, show the some temp result
    % Exit immediately if no corner points were found
    if isempty(scores)
        board = struct('BoardIdx', zeros(3), 'BoardCoords', zeros(3,3,3), ...
                'Energy', Inf, 'isValid', 0);      
        return;
    end

    seedIdx = 1:size(points, 1);
    if size(points, 1) > 100
        % only use corners with high scores as seeds to reduce computation
        if mean(scores) < max(scores) % need this check in case all scores are the same
            seedIdx = seedIdx(scores >= mean(scores));
        end
        [~, sortedIdx] = sort(scores(seedIdx), 'descend');
        seedIdx = seedIdx(sortedIdx);
    end
    
    % grow checkerboard
    board = initializeAndExpandCheckerboard(seedIdx,points,MainPeakFeature,nBins,I);
end

function board = initializeAndExpandCheckerboard(seedIdx,points,MainPeakFeature,nBins,I)
    %% set up parameters
    nSeedPts = length(seedIdx);
    % kdTree for fast NN search
    kdTree = KDTreeSearcher(points);
    % find edge direction
    delta = 2 * pi / nBins;
    epson = cos(deg2rad(30)); % we use larger epson for large distortion
    Angles = 0 + delta .* MainPeakFeature;
    nDirection = 2 * length(Angles);% since DirectVector has 4 directions, but we need 8 neibors
    DirectVector = zeros(nDirection,2);
    for i = 1 : 2 : nDirection
        iAngle = (i + 1) / 2;
        DirectVector(i,:) = [cos(Angles(iAngle)),sin(Angles(iAngle))];
    end
    %% intial and Expand from every seed point
    AllPossibleBoards = cell(nSeedPts,1);
    OptimalNormaliseEnergy = inf;
    OptimalBoardID = 0;
    for i = 1 : nSeedPts
        % search 8 nearest neibors which align the directVector
        [isOk,board2x2] = initial2x2Board(seedIdx(i),DirectVector,kdTree,I);
        if isOk == false
            continue;
        end
        
        % expand the board
        CurrentBoard = board2x2;
        while true && CurrentBoard.Energy < -7.7
            
%             % debug : show current board
%             figure(102);imshow(I);hold on;plot(points(:,1),points(:,2),'ro');hold on;
%             BoardCoords = CurrentBoard.BoardCoords;
%             plot(BoardCoords(:,:,1),BoardCoords(:,:,2),'b*');
            
            [Tmpboard,isStop] = expandBoard(CurrentBoard,DirectVector,kdTree,epson,I);
            if isStop == true
                break;
            else
                CurrentBoard = Tmpboard;
            end
        end
        AllPossibleBoards{i} = CurrentBoard;
%         m = size(CurrentBoard.BoardCoords,1);
%         n = size(CurrentBoard.BoardCoords,2);
        NormalisedEnergy = CurrentBoard.Energy;% / m / n;
        if OptimalNormaliseEnergy > NormalisedEnergy
            OptimalBoardID = i;
            OptimalNormaliseEnergy = NormalisedEnergy;
        end
    end
    %% select the most greatest checkerboard, since we only need 
    if OptimalBoardID > 0
        board = AllPossibleBoards{OptimalBoardID};
    else
        board = struct('BoardIdx', zeros(3,3), 'BoardCoords', zeros(3,3,2),'Energy',...
                    Inf,'MaxTriple',Inf, 'isValid', 0);
    end
end

function [board,isStop] = expandBoard(CurrentBoard,DirectVector,kdTree,epson,I)
% this function expand one step along the boundry
    isStop = false;
    points = kdTree.X;
    
    % expand 4 boundrys
    NewBoundryCoords = cell(4,1);
    NewBoundryIdx = cell(4,1);
    deltaE = zeros(4,1);
    deltaEMaxTriple = zeros(4,1);
    for i = 1 : 4
        
        DirID = 2 * i - 1;
        % expand bundry and compute energy decreasing of this move, from 2 order 
        [deltaE_1,deltaEMaxTriple_1,ThisBoundryMaxTriple_1,NewBoundryCoords_1,NewBoundryIdx_1] = ...
            expandBoundryAndUpdateEnergy(CurrentBoard,DirectVector,DirID,kdTree,epson,I,1);
        [deltaE_2,deltaEMaxTriple_2,ThisBoundryMaxTriple_2,NewBoundryCoords_2,NewBoundryIdx_2] = ...
            expandBoundryAndUpdateEnergy(CurrentBoard,DirectVector,DirID,kdTree,epson,I,2);
        
        % get result from the order which minimize the energy
        if min(min(abs(NewBoundryIdx_1 - NewBoundryIdx_2))) > 0
            deltaE(i) = inf;
            deltaEMaxTriple(i) = inf;
            NewBoundryCoords{i} = [];
            NewBoundryIdx{i} = [];
        else
            if ThisBoundryMaxTriple_1 < ThisBoundryMaxTriple_2
                deltaE(i) = deltaE_1;
                deltaEMaxTriple(i) = deltaEMaxTriple_1;
                NewBoundryCoords{i} = NewBoundryCoords_1;
                NewBoundryIdx{i} = NewBoundryIdx_1;
            else
                deltaE(i) = deltaE_2;
                deltaEMaxTriple(i) = deltaEMaxTriple_2;
                NewBoundryCoords{i} = NewBoundryCoords_2;
                NewBoundryIdx{i} = NewBoundryIdx_2;
            end
        end
    end
    
    % find the boundry which dcrease the energy most
    [mindeltaE,minIdx] = min(deltaE);
    if mindeltaE > 0
        isStop = true;
        board = CurrentBoard;
    else
        % get current board information
        CurBoardIdx = CurrentBoard.BoardIdx;
        CurTotalEnergy = CurrentBoard.Energy;
        [m,n] = size(CurBoardIdx);
        
        % add new boundry to current board
        if minIdx == 1 % direction 1 add to bottom
            BoardIdx = [CurBoardIdx;reshape(NewBoundryIdx{minIdx},1,length(NewBoundryIdx{minIdx}))];
            BoardCoords = reshape(points(BoardIdx,:),m + 1,n,2);
        elseif minIdx == 2 % direction 3 add to left
            BoardIdx = [reshape(NewBoundryIdx{minIdx},length(NewBoundryIdx{minIdx}),1),CurBoardIdx];
            BoardCoords = reshape(points(BoardIdx,:),m,n + 1,2);
        elseif minIdx == 3 % direction 5 add to up
            BoardIdx = [reshape(NewBoundryIdx{minIdx},1,length(NewBoundryIdx{minIdx}));CurBoardIdx];
            BoardCoords = reshape(points(BoardIdx,:),m + 1,n,2);
        else% direction 7 add to right
            BoardIdx = [CurBoardIdx,reshape(NewBoundryIdx{minIdx},length(NewBoundryIdx{minIdx}),1)];
            BoardCoords = reshape(points(BoardIdx,:),m,n + 1,2);
        end
        
        TotalEnergy = CurTotalEnergy + deltaE(minIdx);
        E_maxTriple = deltaEMaxTriple(minIdx);  
        
        % generate new board
        board = struct('BoardIdx', BoardIdx, 'BoardCoords', BoardCoords,'Energy',...
                    TotalEnergy,'MaxTriple',E_maxTriple, 'isValid', 1); 
    end
end

function [deltaE,deltaEMaxTriple,ThisBoundryMaxTriple,NewBoundryCoords,NewBoundryIdx] = ...
    expandBoundryAndUpdateEnergy(Board,DirectVector,DirectID,kdTree,epson,I,Time)

    % initialise
    m = size(Board.BoardCoords,1);
    n = size(Board.BoardCoords,2);
    deltaE = inf;
    deltaEMaxTriple = inf;
    ThisBoundryMaxTriple = inf;
    
    % get boundry and sub-boundry with respect to DirectID 
    % note that this must align the DirectionVector, i.e., Place matrix
    % cjs are exactly boundry coords, cis are rows/cols which close to
    % boundry, and cks are close to cis
    if DirectID == 1
        cis = reshape(Board.BoardCoords(m - 1,:,:),n,2);
        cjs = reshape(Board.BoardCoords(m,:,:),n,2);
        cks = reshape(Board.BoardCoords(m - 2,:,:),n,2);
    elseif DirectID == 3
        cis = reshape(Board.BoardCoords(:,2,:),m,2);
        cjs = reshape(Board.BoardCoords(:,1,:),m,2);
        cks = reshape(Board.BoardCoords(:,3,:),m,2);
    elseif DirectID == 5
        cis = reshape(Board.BoardCoords(2,:,:),n,2);
        cjs = reshape(Board.BoardCoords(1,:,:),n,2);
        cks = reshape(Board.BoardCoords(3,:,:),n,2);
    else
        cis = reshape(Board.BoardCoords(:,n - 1,:),m,2);
        cjs = reshape(Board.BoardCoords(:,n,:),m,2);
        cks = reshape(Board.BoardCoords(:,n - 2,:),m,2);
    end
    
    % search nearist neibors which are algin to DirectVector(DirID) of BoundryCoords
    nBoundrySearchNeibors = 55;
    nBoundryPts = size(cjs,1);
    NewBoundryIdx = zeros(nBoundryPts,1);
    NewBoundryCoords = zeros(nBoundryPts,2);
    BoundryTriples = inf .* ones(nBoundryPts,1);
    BoundryCost = inf .* ones(nBoundryPts,1);
    points = kdTree.X;
    
    % loop for each rest boundry points 
    if Time == 1
        iStart = 1;
        iEnd = nBoundryPts;
        iStep = 1;
    else
        iStart = nBoundryPts;
        iEnd = 1;
        iStep = -1;
    end
    for i = iStart : iStep :iEnd
        CurrentCoord = cjs(i,:);
        SearchRes = knnsearch(kdTree,CurrentCoord,'IncludeTies',true,'K',nBoundrySearchNeibors);
        NIdx = SearchRes{1};
        isSuccesfulExpand = false;
        % find neibor which align the DirectVector
        for j = 2 : length(NIdx)
            thisID = NIdx(j);
            thisCoord = points(thisID,:);
            curToThisVector = thisCoord - CurrentCoord;
            
            % regular
            isRegular = false;
            if abs(norm(cks(i,:) - cis(i,:)) - norm(cjs(i,:) - cis(i,:))) < 20
%                 maxGridSize = max(norm(cks(i,:) - cis(i,:)),norm(cjs(i,:) - cis(i,:),norm(curToThisVector)));
%                 minGridSize = min(norm(cks(i,:) - cis(i,:)),norm(cjs(i,:) - cis(i,:),norm(curToThisVector)));
%                 if maxGridSize - minGridSize < 20       
%                 end
                isRegular = true;
            elseif norm(cks(i,:) - cis(i,:)) > norm(cjs(i,:) - cis(i,:))
                if norm(cjs(i,:) - cis(i,:)) > norm(curToThisVector)
                    isRegular = true;
                end
            elseif norm(cks(i,:) - cis(i,:)) < norm(cjs(i,:) - cis(i,:))
                if norm(cjs(i,:) - cis(i,:)) < norm(curToThisVector)
                    isRegular = true;
                end
            end
            
            if isRegular == false
                continue;
            end
            
            curToThisVector = curToThisVector ./ norm(curToThisVector);
            dotDiff = dot(curToThisVector, DirectVector(DirectID,:));
            if dotDiff > epson
                score = I(thisCoord(2),thisCoord(1));% here we can change
                FlatTriple = norm(cis(i,:) + thisCoord - 2 * cjs(i,:)) / norm(cis(i,:) - thisCoord);
                if i == iStart
                    crossId = i + iStep * 2;
                    midId = i + iStep;
                    Cross = norm(cis(crossId,:) + thisCoord - 2 * cjs(midId ,:)) / norm(cis(crossId,:) - thisCoord);
                    ThisCost = FlatTriple + Cross;
                else
                    e1 = thisCoord - NewBoundryCoords(i - iStep,:);
                    e2 = cjs(i,:) - cjs(i - iStep,:);
                    e3 = NewBoundryCoords(i - iStep,:) - cjs(i - iStep,:);
                    e4 = thisCoord - cjs(i,:);
                    e5 = thisCoord - cjs(i - iStep,:);
                    ThisCost = norm(e1 + e4 - e5) ./ norm(e1 - e4) + norm(e2 + e3 - e5) ./ norm(e2 - e3)...
                            + FlatTriple;
                end
                ThisCost = ThisCost / score;
                
                % among all neibors, we find the one which minimise the triple
                isNotRepeat = true;
                for id = 1 : nBoundryPts
                    if thisCoord == NewBoundryCoords(id,:)
                        isNotRepeat = false;
                    end
                end
                % find the best triple which is not contain repeat
                if BoundryCost(i) > ThisCost && isNotRepeat == true
                    NewBoundryIdx(i) = thisID;
                    NewBoundryCoords(i,:) = thisCoord;
                    BoundryTriples(i) = FlatTriple;
                    BoundryCost(i) = ThisCost;
                    isSuccesfulExpand = true;
                end
            end
            
        end
        % count how many boundry points are successfully expand
        if isSuccesfulExpand == false
            return;
        end
    end
    
    % compute Energy of this move
    % compute deltaEMaxTriple
    deltaEMaxTriple = Board.MaxTriple;
    ThisBoundryMaxTriple = max(BoundryTriples);
    for i = 1 : nBoundryPts - 2
        tmpTriple = norm(NewBoundryCoords(i,:) + NewBoundryCoords(i + 2,:) - 2 * NewBoundryCoords(i + 1,:))...
                    ./norm(NewBoundryCoords(i,:) - NewBoundryCoords(i + 2,:));
        ThisBoundryMaxTriple = max(ThisBoundryMaxTriple,tmpTriple);
    end
    % compute deltaEStruct
    deltaEMaxTriple = max(ThisBoundryMaxTriple,deltaEMaxTriple);
    deltaEStruct = (m * n + nBoundryPts) * deltaEMaxTriple - (m * n) * Board.MaxTriple;
    deltaEcorner = -1 * nBoundryPts;
    deltaE = deltaEcorner + deltaEStruct;
    % deltaEStruct can not be too much
    if abs(deltaE) < deltaEStruct
        deltaE = inf;
    end
end


function [isOk,board] = initial2x2Board(seedIdx,DirectVector,kdTree,I) 
    % setup
    isOk = true;
    nSearchNeibors = 50;
    epson = cos(deg2rad(25));% strict than boundry expand
    points = kdTree.X; 
    nDirection = size(DirectVector,1);
    SeedCoord = points(seedIdx,:);
    SeedNeiborsIdx = zeros(nDirection,1);
    SeedNeiborsCoords = zeros(nDirection,2);
    
%     % debug : show vector
%     figure(100);imshow(I);hold on;plot(points(:,1),points(:,2),'ro');hold on;
%     seedPts = points(seedIdx,:);
%     ePts(:,1) = seedPts(1) + 50 .* DirectVector(:,1);
%     ePts(:,2) = seedPts(2) + 50 .* DirectVector(:,2);
%     for ii = 1 : 2 : size(DirectVector,1)
%         plot([seedPts(1),ePts(ii,1)],[seedPts(2),ePts(ii,2)],'r-');
%     end
%     for ii = 1 : 2 : size(DirectVector,1)
%         label = sprintf('%d',ii);
%         text(ePts(ii,1), ePts(ii,2), label, 'BackgroundColor', [1 1 1]);
%     end
    
    % search nearist neibor
    SearchRes = knnsearch(kdTree,SeedCoord,'IncludeTies',true,'K',nSearchNeibors);
    NeiborIdx = SearchRes{1};
    % set all direction are avaliable
    isAvaliable = zeros(nDirection,1);
    isAvaliable([1,3,5,7]) = 1;
    Cost = inf * ones(nDirection,1);
    
    % find neirbors which align the directVector, 4 neibors
    for j = 2 : length(NeiborIdx)
        thisID = NeiborIdx(j);
        thisCoord = points(thisID,:);
        seedToThisVector = thisCoord - SeedCoord;
        seedToThisVector = seedToThisVector ./ norm(seedToThisVector);
        % test if this neibor is align with the directVector, only for 1,3,5,7
        for iDir = 1 : 2 : nDirection
            dotDiff = dot(seedToThisVector, DirectVector(iDir,:));
            if dotDiff > epson
                thisCost = norm(thisCoord - SeedCoord) / I(thisCoord(2),thisCoord(1));
                % make sure there do not contain any duplex element
                isNotRepeat = true;
                for id = 1 : size(SeedNeiborsCoords,1)
                    if thisCoord == SeedNeiborsCoords(id,:)
                        isNotRepeat = false;
                    end
                end
                if Cost(iDir) > thisCost && isNotRepeat == true
                    SeedNeiborsIdx(iDir,:) = thisID;
                    SeedNeiborsCoords(iDir,:) = thisCoord;
                    Cost(iDir) = thisCost;
                    isAvaliable(iDir) = 0;
                end
            end
        end
    end

%     % only for debug
%     figure(101);imshow(I);hold on;plot(points(:,1),points(:,2),'ro');hold on;
%     plot(SeedCoord(1),SeedCoord(2),'g+');hold on;
%     plot(SeedNeiborsCoords(1:2:8,1),SeedNeiborsCoords(1:2:8,2),'b*');
    
    % if can not found 4 correct neibors of seed, then remove this seed point
    if sum(isAvaliable) > 0
        isOk = false;      
        board = struct('BoardIdx', zeros(3,3), 'BoardCoords', zeros(3,3,2),'Energy',...
                    Inf,'MaxTriple',Inf, 'isValid', 0); 
        return;
    end
    
    % find another 4 neibors for 45 direction
    k = 4;
    for iDir = 2 : 2 : nDirection
        minDist = inf;
        PrevCoord = SeedNeiborsCoords(iDir - 1,:);
        if iDir >= 8
            NextID = 1;
        else
            NextID = iDir + 1;
        end
        NextCoord = SeedNeiborsCoords(NextID,:);
        % knn search for neibors
        SearchRes = knnsearch(kdTree,PrevCoord,'IncludeTies',true,'K',nSearchNeibors);
        IDs = SearchRes{1};
        for j = 2 : length(IDs)
            thisID = IDs(j);
            thisCoord = points(thisID,:);
            thisToPrevVector = thisCoord - PrevCoord;
            thisToPrevVector = thisToPrevVector ./ norm(thisToPrevVector);
            % test if this neibor is align with the CurrentCoord, 
            if dot(thisToPrevVector, DirectVector(NextID,:)) > epson
                e1 = PrevCoord - SeedCoord;
                e2 = NextCoord - SeedCoord;
                e3 = thisCoord - SeedCoord;
                thisDist = norm(e1 + e2 - e3) / norm(e1 - e2) / I(thisCoord(2),thisCoord(1));%
                % make sure there do not contain any duplex element
                isNotRepeat = true;
                for id = 1 : size(SeedNeiborsCoords,1)
                    if thisCoord == SeedNeiborsCoords(id,:)
                        isNotRepeat = false;
                    end
                end
                if thisDist < minDist && isNotRepeat == true
                    SeedNeiborsIdx(iDir,:) = thisID;
                    SeedNeiborsCoords(iDir,:) = thisCoord;
                    k = k + 1;
                    minDist = thisDist;
                end
            end
        end
    end
    % if can not found 8 correct neibors of seed, then remove this seed point
    if k < nDirection
        isOk = false;      
        board = struct('BoardIdx', zeros(3,3), 'BoardCoords', zeros(3,3,2),'Energy',...
                    Inf,'MaxTriple',Inf, 'isValid', 0); 
        return;
    end
    

    % compute energy of initial 2x2 board
    AllPointsIn2x2Board = [SeedNeiborsCoords;SeedCoord];
    AllPointsIn2x2Board_Idx = [SeedNeiborsIdx;seedIdx];
    [board,E] = compute2x2BoardEnergy(AllPointsIn2x2Board,AllPointsIn2x2Board_Idx);
    
%     % only for debug
%     figure(101);imshow(I);hold on;plot(points(:,1),points(:,2),'ro');hold on;
%     plot(SeedCoord(1),SeedCoord(2),'g+');hold on;
%     plot(SeedNeiborsCoords(:,1),SeedNeiborsCoords(:,2),'b*');

    if E >= -7.5
        isOk = false;
        return;% if initial board is terrible then switch to next seed point
    end
end

function [board,E] = compute2x2BoardEnergy(AllPointsIn2x2Board,AllPointsIn2x2Board_Idx)

    Place = [4,5,6;
             3,9,7;
             2,1,8];

    nDirection = 8;
    E_corners = -1 * (1 + nDirection);
    E_maxTriple = 0;
    for iRow = 1 : 3
        ci = AllPointsIn2x2Board(Place(iRow,1),:);
        cj = AllPointsIn2x2Board(Place(iRow,2),:);
        ck = AllPointsIn2x2Board(Place(iRow,3),:);
        E_maxTriple = max(E_maxTriple,norm(ci + ck - 2 * cj) / norm(ci - ck));
    end
    for iCol = 1 : 3
        ci = AllPointsIn2x2Board(Place(1,iCol),:);
        cj = AllPointsIn2x2Board(Place(2,iCol),:);
        ck = AllPointsIn2x2Board(Place(3,iCol),:);
        E_maxTriple = max(E_maxTriple,norm(ci + ck - 2 * cj) / norm(ci - ck));
    end
    E_struct = E_maxTriple * (1 + nDirection);
    E = E_corners + E_struct;
    % 2x2 board structure
    BoardIdx = AllPointsIn2x2Board_Idx(Place);
    BoardCoords = reshape(AllPointsIn2x2Board(Place,:),3,3,2);
    board = struct('BoardIdx', BoardIdx, 'BoardCoords', BoardCoords,'Energy',E,...
                    'MaxTriple',E_maxTriple, 'isValid', true); 
end