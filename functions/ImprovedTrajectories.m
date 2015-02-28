function F = ImprovedTrajectories
    F.RunImprovedTrajectories = @RunImprovedTrajectories;
    F.AnalyzeOutput = @AnalyzeOutput;
    F.GroupTrajectoriesManual = @GroupTrajectoriesManual;
    F.GroupTrajectoriesKmeans = @GroupTrajectoriesKmeans;
    F.GetSeedTrjs = @GetSeedTrjs;
    F.PlotSeeds = @PlotSeeds;
    F.PlotGroups = @PlotGroups;
    F.PlotAllTrajectories = @PlotAllTrajectories;    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Public Functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [status,exeOutput] = RunImprovedTrajectories(videoPath)        
    [status,exeOutput] = system(['"C:\Users\ckocb_000\Documents\Visual Studio 2012\Projects\ImprovedTrajectories\x64\Debug\ImprovedTrajectories.exe" "' videoPath '"']);
end

function trajectories = AnalyzeOutput(exeOutput)

    if isempty(exeOutput); error('ImpTrajectories:propertyCheck', 'First run ''RunImpTrajectories''.'); end

    tempTrajectories = regexp(exeOutput,'\n','split');
    numberOfTrajectories = size(tempTrajectories,2);

    Util = UtilFunctions;
    for k = 1:numberOfTrajectories-2 % the last element is empty, and the previous one contains info about the movie

        currentTempTrajectory = regexp(tempTrajectories{k},'\s','split');
        currentTempTrajectory = Util.FlattenCellArray(currentTempTrajectory);
        currentTempTrajectory = str2double(currentTempTrajectory);     

        currentTrajectory = struct();

        currentTrajectory.frameNum  = currentTempTrajectory(1);
        currentTrajectory.mean_x    = currentTempTrajectory(2);   
        currentTrajectory.mean_y    = currentTempTrajectory(3);    
        currentTrajectory.var_x     = currentTempTrajectory(4); 
        currentTrajectory.var_y     = currentTempTrajectory(5); 
        currentTrajectory.length    = currentTempTrajectory(6); 
        currentTrajectory.scale     = currentTempTrajectory(7);     
        currentTrajectory.x_pos     = currentTempTrajectory(8);       
        currentTrajectory.y_pos     = currentTempTrajectory(9);       
        currentTrajectory.t_pos     = currentTempTrajectory(10);  

        trajectoryLenght = floor((size( currentTempTrajectory, 2 ) - 10)/2);%-96 -108 -96 -96; 
        if trajectoryLenght > 1
%             currentTrajectory.trajectory=  currentTempTrajectory( 11 : 10+trajectoryLenght ); 
            currentTrajectory.trajectory= reshape(  currentTempTrajectory( 11 : 10+trajectoryLenght*2 ) , [2 trajectoryLenght] ); 
        else
            currentTrajectory.trajectory = [];
        end
        
%         currentTrajectory.HOG       =  currentTempTrajectory( trajectoryLenght+11  : trajectoryLenght+106 ) ; %(default 96 dimension)        
%         currentTrajectory.HOF       =  currentTempTrajectory( trajectoryLenght+107 : trajectoryLenght+214 ) ; %(default 108 dimension) 
%         currentTrajectory.MBHx      =  currentTempTrajectory( trajectoryLenght+215 : trajectoryLenght+310 ) ; %(default 96 dimension)
%         currentTrajectory.MBHy      =  currentTempTrajectory( trajectoryLenght+311 : trajectoryLenght+406 ) ; %(default 96 dimension)
        
%         currentTrajectory.HOG       = reshape(  currentTempTrajectory( trajectoryLenght+11  : trajectoryLenght+106 ) , [8 2 2 3] ); %(default 96 dimension)        
%         currentTrajectory.HOF       = reshape(  currentTempTrajectory( trajectoryLenght+107 : trajectoryLenght+214 ) , [9 2 2 3] ); %(default 108 dimension) 
%         currentTrajectory.MBHx      = reshape(  currentTempTrajectory( trajectoryLenght+215 : trajectoryLenght+310 ) , [8 2 2 3] ); %(default 96 dimension)
%         currentTrajectory.MBHy      = reshape(  currentTempTrajectory( trajectoryLenght+311 : trajectoryLenght+406 ) , [8 2 2 3] ); %(default 96 dimension)

        trajectories(k) = currentTrajectory;
    end

end

function seeds = GetSeedTrjs(trajectories,keyFrames,shotBoundaries,staticSaliency,adjacencyDistance)

    if nargin < 5
        adjacencyDistance = 10;
    end

    nFrames = trajectories(end).frameNum;
    nShots = size(shotBoundaries,1)-1;
    nTrajectories = size(trajectories,2);
    
    staticSaliency(staticSaliency<0.7) = 0;
  
    % Group trajectories according to their start frame
    trajectoryLengths = arrayfun(@(X)size(trajectories(X).trajectory,2),1:nTrajectories);
    trajectoryEnds = [trajectories.frameNum];
    trajectoryStarts = trajectoryEnds - trajectoryLengths + 1;
    
    % Get valid trajectories (remove inter-shot trajectories)
    isInsideShotBoundaries = arrayfun(@(X)(trajectoryStarts > shotBoundaries(X) & ...
        trajectoryEnds < shotBoundaries(X+1)),1:nShots,'UniformOutput',false);
    isInsideShotBoundaries = any(cell2mat(isInsideShotBoundaries'),1);
    validTrajectoryIndices = find(isInsideShotBoundaries)';
    trajectoryStarts = trajectoryStarts(isInsideShotBoundaries)'; 
    nTrajectories = size(trajectories,2);
    
    % Get seed trajectories from key frames
    seeds = [];
    seedStarts = [];
    for k = 1:size(keyFrames,1)
        
        currentFrame = keyFrames(k);
        
        % Get trajectories that can touch current key frame
        possibleTrajs = find(trajectoryStarts>=currentFrame-12 & ...
            trajectoryStarts<=currentFrame);
        possibleStarts = trajectoryStarts(possibleTrajs);
        possibleTrajs = validTrajectoryIndices(possibleTrajs);
        nrIndices = size(possibleTrajs,1);
        
        % Get positions of these trajectories on current key frame
        currentX = arrayfun(@(X)trajectories(possibleTrajs(X)).trajectory( ...
            2,currentFrame-possibleStarts(X)+1),1:nrIndices);
        currentY = arrayfun(@(X)trajectories(possibleTrajs(X)).trajectory( ...
            1,currentFrame-possibleStarts(X)+1),1:nrIndices);
        currentX = floor(currentX)';
        currentY = floor(currentY)';
        currentX(currentX==0) = 1;
        currentY(currentY==0) = 1;
        
        % Select trajectories that pass through salient regions on keyframe
        salientTrajectories = cell2mat(arrayfun(@(X)staticSaliency(currentY(X),currentX(X),k),...
            1:nrIndices,'UniformOutput',false));
        seeds = [seeds;possibleTrajs(find(salientTrajectories))];
        seedStarts = [seedStarts;possibleStarts(find(salientTrajectories))];
        
    end

    % Find seed trajectories from ramaining frames
    for t = 1:2
        for k = 1:nFrames

            % Get trajectories starting at current frame
            possibleTrajs = validTrajectoryIndices(trajectoryStarts==k);
            possibleTrajsStarts = trajectoryStarts(trajectoryStarts==k);
            nrTrajs = size(possibleTrajs,1);
            if ~nrTrajs; continue; end;

            % Get trajectory data of found trajectories
            possibleTrajsData   = [trajectories(possibleTrajs).trajectory]';
            possibleTrajsX      = reshape(possibleTrajsData(:,1),[13,size(possibleTrajsData,1)/13])'; 
            possibleTrajsY      = reshape(possibleTrajsData(:,2),[13,size(possibleTrajsData,1)/13])'; 

            % Get seeds from upcoming 10 frames
            possibleSeeds       = arrayfun(@(X)seeds(seedStarts==(k+X)),0:9, 'UniformOutput', false);
            possibleSeedStarts  = arrayfun(@(X)seedStarts(seedStarts==(k+X)),0:9, 'UniformOutput', false);
            possibleSeeds       = cell2mat(possibleSeeds(:));
            possibleSeedStarts  = cell2mat(possibleSeedStarts(:));
            nrSeeds             = size(possibleSeeds,1);
            if ~nrSeeds; continue; end;

            % Get trajectory data of seeds
            possibleSeedsData   = [trajectories(possibleSeeds).trajectory]';
            possibleSeedsX      = reshape(possibleSeedsData(:,1),[13,size(possibleSeedsData,1)/13])'; 
            possibleSeedsY      = reshape(possibleSeedsData(:,2),[13,size(possibleSeedsData,1)/13])';

            % Get the difference between seed starts and current frame
            possibleStartDiff = possibleSeedStarts-k;

            % Get overlaping portions of seeds with trajectories
            possibleSeedsX = cell2mat(arrayfun(@(X)[ possibleSeedsX(X,1:end-possibleStartDiff(X))'; ...
                zeros(possibleStartDiff(X),1)], 1:nrSeeds,'uni',false))';
            possibleSeedsY = cell2mat(arrayfun(@(X)[ possibleSeedsY(X,1:end-possibleStartDiff(X))'; ...
                zeros(possibleStartDiff(X),1)], 1:nrSeeds,'uni',false))';

            possibleSeedsX = repmat(possibleSeedsX,[1,1,nrTrajs]);
            possibleSeedsY = repmat(possibleSeedsY,[1,1,nrTrajs]);

            possibleSeedsX = permute(possibleSeedsX,[3 2 1]);
            possibleSeedsY = permute(possibleSeedsY,[3 2 1]);

            % Calculate the distance of current trajectories with seeds
            possibleTrajsX = repmat(possibleTrajsX,[1,1,nrSeeds]);
            possibleTrajsY = repmat(possibleTrajsY,[1,1,nrSeeds]);

            possibleTrajsX = cell2mat(arrayfun(@(X)[possibleTrajsX(:,possibleStartDiff(X)+1:end,1),...
                zeros(nrTrajs,possibleStartDiff(X))], 1:nrSeeds,'uni',false)');
            possibleTrajsY = cell2mat(arrayfun(@(X)[possibleTrajsY(:,possibleStartDiff(X)+1:end,1),...
                zeros(nrTrajs,possibleStartDiff(X))], 1:nrSeeds,'uni',false)');

            possibleTrajsX = reshape(possibleTrajsX,[nrTrajs nrSeeds 13]);
            possibleTrajsY = reshape(possibleTrajsY,[nrTrajs nrSeeds 13]);

            possibleTrajsX = permute(possibleTrajsX,[1 3 2]);
            possibleTrajsY = permute(possibleTrajsY,[1 3 2]);

            distanceX = possibleTrajsX - possibleSeedsX;
            distanceY = possibleTrajsY - possibleSeedsY;

            distanceX(distanceX==0) = nan;
            distanceX = nanmean(distanceX,2);
            distanceY(distanceY==0) = nan;
            distanceY = nanmean(distanceY,2);

            distance = sqrt(distanceX.^2+distanceY.^2); 
            
            [trajInd,~] = ind2sub([nrTrajs nrSeeds], find(distance < adjacencyDistance));

            seeds = [seeds;possibleTrajs(unique(trajInd,'rows'))];
            seedStarts = [seedStarts;possibleTrajsStarts(unique(trajInd,'rows'))];

        end 
    end
   
    seeds = unique(seeds,'rows');

end

function groups = GroupTrajectoriesManual(trajectories, shotBoundaries, adjacencyDistance)

    % Group trajectories according to their starting frame
    nFrames = trajectories{end}.frameNum;
    accordingToStart = cell(nFrames,1);
    for i = 1:nFrames; accordingToStart{i} = []; end
    for k = 1:size(trajectories,1)

        currentTrajectory = trajectories{k};

        trajectoryLength = size(currentTrajectory.trajectory,2);
        trajectoryEnd = currentTrajectory.frameNum;
        trajectoryStart = trajectoryEnd - trajectoryLength + 1;

        if trajectoryStart > shotBoundaries(1) && trajectoryEnd < shotBoundaries(2)
            accordingToStart{trajectoryStart - shotBoundaries(1)}(end+1) = k;
        end

    end

    % Find adjacent trajectories
    adjacencyDistance = 10;
    adjacencies = zeros(100000,2);
    adjacencies(1,1) = 2; % keep the last index in the first element
    for k = 1:nFrames

        currentStart = accordingToStart{k};
        for i = 1:10:size(currentStart,2)

            traji = trajectories{currentStart(i)}.trajectory; 

            % Check trajectories starting from the same frame
            for j = i+1:size(currentStart,2)
                trajj = trajectories{currentStart(j)}.trajectory;
                distance = abs(mean(traji(2,:) - trajj(2,:)));
                if distance < adjacencyDistance
                    adjacencies(adjacencies(1,1),:) = [currentStart(i),currentStart(j)];
                    adjacencies(1,1) = adjacencies(1,1) + 1;
                end
            end

            % Check trajectories starting from the -1 frame
            if k > 5
                for j = 1:5
                    possibleStart = accordingToStart{k-j};
                    for t = 1:size(possibleStart,2)
                        trajj = trajectories{possibleStart(t)}.trajectory;
                        distancex = mean(abs(traji(2,1:end-j) - trajj(2,j+1:end)));
                        distancey = mean(abs(traji(1,1:end-j) - trajj(1,j+1:end)));
                        distance = sqrt(distancex^2+distancey^2);
                        if distance < adjacencyDistance
                            adjacencies(adjacencies(1,1),:) = [currentStart(i),possibleStart(t)];
                            adjacencies(1,1) = adjacencies(1,1) + 1;
                        end
                    end
                end
            end
        end
    end 
    adjacencies(~any(adjacencies,2),:) = []; % remove zero rows

    % Group adjacent trajectories
    GROUP_NUMBER = 1000;
    groups = zeros(2000,GROUP_NUMBER);
    groupIndex = 1;
    groupItemIndex = ones(GROUP_NUMBER,1);
    while 1

        if size(adjacencies,1) == 1; break; end

        % Find all adjacent trajectories with the current one
        currentTrajectoryIndex = adjacencies(end,1);
        indices = find((adjacencies == currentTrajectoryIndex));
        [rowind,colind]=ind2sub(size(adjacencies),indices);
        adjacentIndices = sub2ind(size(adjacencies),rowind,mod((colind),2)+1);
        rowind = [rowind;rowind(1)];
        colind = [colind;colind(1)];
        adjacentIndices = [adjacentIndices;indices(1)];

        % Check if any of the adjacent trajectory is already member of a group
        currentGroup = 0;
        for i = 1:size(adjacentIndices)
            currentGroup = find(groups == adjacencies(rowind(i),colind(i)));
            if ~isempty(currentGroup)
                [~,currentGroup] = ind2sub(size(groups),currentGroup(1));
                break;
            end
        end

        % If non of the adjacent trajectories is a part of a group, create a
        % new group for them
        if isempty(currentGroup)
            currentGroup = groupIndex;
            groupIndex = groupIndex +1;
        end

        % Assign trajectories to the defined group
        currentGroupItemIndex = groupItemIndex(currentGroup);
        groups(currentGroupItemIndex:currentGroupItemIndex +  ...
            size(adjacentIndices,1)-1,currentGroup) = adjacencies(adjacentIndices);
        groupItemIndex(currentGroup) = currentGroupItemIndex+size(adjacentIndices,1);

        % Remove the processed trajectories
        adjacencies(rowind,:)=[];

    end
    groups(:,~any(groups,1)) = [];  % Remove all zero columns
    groups(~any(groups,2),:) = [];  % Remove all zero rows
end

function groups = GroupTrajectoriesKmeans(trajectories, shotBoundaries)
    
    forKMeans = zeros(size(trajectories,1),6); %[meanX,meanY,firstD,secondD]
    for k = 1:size(trajectories,1)

        currentTrajectory = trajectories{k};

        trajectoryLength = size(currentTrajectory.trajectory,2);
        trajectoryEnd = currentTrajectory.frameNum;
        trajectoryStart = trajectoryEnd - trajectoryLength + 1;

        if trajectoryStart > shotBoundaries(1) && trajectoryEnd < shotBoundaries(2)   
            tr = currentTrajectory.trajectory;
    %         forKMeans(k,:) = dtw_c(tr(1,:),tr(2,:));
    %         X = fft(tr(1,:));
    %         Y = fft(tr(2,:));
    %         forKMeans(kMeansIndex,:) = [real(X(1)),real(Y(1))];  
    %         correspondance(kMeansIndex,:) = [kMeansIndex,k];

            currentTrajectory.distance = sqrt(tr(1,:).^2 + tr(2,:).^2);
            meanDistance = mean(currentTrajectory.distance);
            p = polyfit(1:trajectoryLength,currentTrajectory.distance,3);
            pd = polyder(p);
            pdd = polyder(pd);

            forKMeans(k,:) = [meanDistance,pd,pdd]; 
        else
            forKMeans(k,:) = [0,0,0,0,0,0];
        end
    end
    forKMeans(~any(forKMeans,2),:) = [];

    GROUP_NUMBER = 5;
    groups = kmeans(forKMeans,GROUP_NUMBER);

end

function PlotSeeds(trajectories,seeds)

    for k = 1:size(seeds,1)

        currentTrajectory = trajectories(seeds(k));
        tr = currentTrajectory.trajectory;
        trajectoryLength = size(currentTrajectory.trajectory,2);
        trajectoryEnd = currentTrajectory.frameNum;
        trajectoryStart = trajectoryEnd - trajectoryLength + 1;

        hold on;
        plot3(tr(2,:),(trajectoryStart-1):(trajectoryEnd-1),...
            tr(1,:),'LineWidth',2);

    end
%     xlim([1 320])
%     ylim([1 114])
%     zlim([1 240])
%     ylim([1 (trajectoryEnd-trajectoryStart)])
%     set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
    xlabel('Video Width','FontWeight','Bold');
    ylabel('Time / Frames','FontWeight','Bold');
    zlabel('Video Height','FontWeight','Bold');
    view(3);
    
end

function PlotGroups(trajectories,groups,shotBoundaries)

    colorArr = {'b','m','k','g','r','c','y'}; 
%     markerArr = {'o','d','+','*','x','.','s','p'}; 

    for k = 1:size(groups,2)

        currentGroup = groups(:,k);
        currentGroup(~any(currentGroup,2),:) = [];
        if size(currentGroup,1) < 300; continue; end;
        
        for i = 1: 50 : size(currentGroup,1)
            
            currentTrajectory = trajectories{currentGroup(i)};
            tr = currentTrajectory.trajectory;

            trajectoryLength = size(tr,2);

            trajectoryEnd = currentTrajectory.frameNum;
            trajectoryStart = trajectoryEnd - trajectoryLength + 1;
            hold on;
            plot3(tr(2,:),(trajectoryStart-1):(trajectoryEnd-1),...
                tr(1,:),'Color',colorArr{mod(k,7)+1},'LineWidth',2);
        end
    end
    set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
    xlabel('Video Width','FontWeight','Bold');
    ylabel('Time / Frames','FontWeight','Bold');
    zlabel('Video Height','FontWeight','Bold');
    axis tight;
    view(3);
    toc;

end

function PlotAllTrajectories(trajectories,shotBoundaries)

    h1 = figure;
    for k = 1:100:size(trajectories,2)

        currentTrajectory = trajectories(k);
        tr = currentTrajectory.trajectory;

        trajectoryLength = size(tr,2);

        trajectoryEnd = currentTrajectory.frameNum;
        trajectoryStart = trajectoryEnd - trajectoryLength + 1;
        if (trajectoryEnd > shotBoundaries(2)) || (trajectoryStart < shotBoundaries(1)); continue; end

        hold on;
%         plot3(tr(2,:),(trajectoryStart-1)*10:10:(trajectoryEnd-1)*10,tr(1,:),'r')
        plot(tr(2,:),(trajectoryStart-1)*10:10:(trajectoryEnd-1)*10,'Color','r','linewidth',3)
    end
    print(h1,'-djpeg','-r150','aerialtrjs.jpg');
end

% function trajectoriesByFrame = GetTrajectoriesByFrame(trajectories)
% 
%     if isempty(trajectories); error('ImpTrajectories:propertyCheck', 'First run ''AnalyzeOutput'' function.'); end
% 
%     numberOfTrajectories = size(trajectories,2);
%     numberOfFrames = trajectories(end).frameNum + 1;
%     trajectoriesByFrame = cell(numberOfFrames,1);
% 
%     for i = 1: numberOfTrajectories
% 
%         currentTrajectory = trajectories(i);
%         trajectoryLength = size(currentTrajectory.trajectory,2);
%         
% %         if trajectoryLength < 25; return; end
%         
%         trajectoryEnd = currentTrajectory.frameNum;
%         trajectoryStart = trajectoryEnd - trajectoryLength + 1; 
% 
%         for k = 1 : trajectoryLength
% 
%             currentFrame = trajectoryStart + k ;
%             descriptor = [ currentTrajectory.trajectory(1,k); ...
%                            currentTrajectory.trajectory(2,k) ];
% 
% %             if ~isempty(trajectoriesByFrames{currentFrame}); trajectoriesByFrames{currentFrame} = []; end
%             trajectoriesByFrame{currentFrame} = [trajectoriesByFrame{currentFrame} , descriptor];
% 
%         end
% 
%     end
% 
% end