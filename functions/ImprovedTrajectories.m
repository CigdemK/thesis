function F = ImprovedTrajectories
    F.RunImprovedTrajectories = @RunImprovedTrajectories;
    F.AnalyzeOutput = @AnalyzeOutput;
    F.GetTrajectoriesByFrame = @GetTrajectoriesByFrame;
    F.GroupTrajectories = @GroupTrajectories;
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

    trajectories = cell(numberOfTrajectories-2,1);
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

        trajectories{k} = currentTrajectory;
    end

end

function trajectoriesByFrame = GetTrajectoriesByFrame(trajectories)

    if isempty(trajectories); error('ImpTrajectories:propertyCheck', 'First run ''AnalyzeOutput'' function.'); end

    numberOfTrajectories = size(trajectories,1);
    numberOfFrames = trajectories{end}.frameNum + 1;
    trajectoriesByFrame = cell(numberOfFrames,1);

    for i = 1: numberOfTrajectories

        currentTrajectory = trajectories{i};
        trajectoryLength = size(currentTrajectory.trajectory,2);
        
        if trajectoryLength < 25; return; end
        
        trajectoryEnd = currentTrajectory.frameNum;
        trajectoryStart = trajectoryEnd - trajectoryLength + 1; 

        for k = 1 : trajectoryLength

            currentFrame = trajectoryStart + k ;
            descriptor = [ currentTrajectory.trajectory(1,k); ...
                           currentTrajectory.trajectory(2,k) ];

%             if ~isempty(trajectoriesByFrames{currentFrame}); trajectoriesByFrames{currentFrame} = []; end
            trajectoriesByFrame{currentFrame} = [trajectoriesByFrame{currentFrame} , descriptor];

        end

    end

end

function groups = GroupTrajectories(trajectories, saliencyMap, shotBoundaries)

    if nargin < 4
        ADJACENCY_DISTANCE = 10;
    end

    nFrames = trajectories{end}.frameNum;
  
    % Group trajectories according to their start frame
    tic;
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
        
        % Find l1 distance of x,y points to the left-top corner
        tr = currentTrajectory.trajectory;
        currentTrajectory.distance = sqrt(tr(1,:).^2 + tr(2,:).^2);

    end
    toc;

    % Find seed trajectories that starts from first frame
    seeds = zeros(10000,nFrames);
    seeds(1,:) = 2;

    currentStart = accordingToStart{1};
    for i = 1:size(currentStart,2)
        traj = trajectories{currentStart(i)}; 
        if saliencyMap(round(traj.mean_y),round(traj.mean_x))
            seeds(seeds(1,1),1) = currentStart(i);
            seeds(1,1) = seeds(1,1) + 1;
        end
    end

    % Find seed trajectories from ramaining frames
    for k = 2:size(accordingToStart,1)
        currentStart = accordingToStart{k};
        for i = 1:size(currentStart,2)

            traj = trajectories{currentStart(i)}.trajectory; 

            % Check seeds of the previous frame
            for t = 2:seeds(1,k-1)-1
                seed = trajectories{seeds(t,k-1)}.trajectory;
                distance = mean(sqrt((traj(2,1:end-1) - seed(2,2:end)).^2 + ...
                    (traj(1,1:end-1) - seed(1,2:end)).^2 ));
                if distance < ADJACENCY_DISTANCE
                    seeds(seeds(1,k),k) = currentStart(i);
                    seeds(1,k) = seeds(1,k) + 1;
                    break;
                end
            end

        end
        toc;
    end
    seeds(1,:) = [];
    seeds = seeds(:);
    seeds(~any(seeds,2),:) = [];
    toc;
    
    %Plot seeds
    for k = 1:5:size(seeds,1)

        currentTrajectory = trajectories{seeds(k)};
        tr = currentTrajectory.trajectory;
        trajectoryLength = size(currentTrajectory.trajectory,2);
        trajectoryEnd = currentTrajectory.frameNum;
        trajectoryStart = trajectoryEnd - trajectoryLength + 1;

        hold on;
        plot3(tr(2,:),(trajectoryStart-1):(trajectoryEnd-1),...
            tr(1,:),'LineWidth',2);

    end
    xlim([1 size(saliencyMap,2)])
%     ylim([1 (trajectoryEnd-trajectoryStart)])
    set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
    xlabel('Video Width','FontWeight','Bold');
    ylabel('Time / Frames','FontWeight','Bold');
    zlabel('Video Height','FontWeight','Bold');
    view(3);
    
%     % Find seed trajectories from ramaining frames
%     for k = 2:size(accordingToStart,1)
%         currentStart = accordingToStart{k};
%         for i = 1:size(currentStart,2)
% 
%             traj = trajectories{currentStart(i)}.trajectory; 
% 
%             % Check seeds of the previous frame
%             for t = 2:seeds(1,k-1)-1
%                 seed = trajectories{seeds(t,k-1)}.trajectory;
%                 distance = mean(sqrt((traj(2,1:end-1) - seed(2,2:end)).^2 + ...
%                     (traj(1,1:end-1) - seed(1,2:end)).^2 ));
%                 if distance < ADJACENCY_DISTANCE
%                     seeds(seeds(1,k),k) = currentStart(i);
%                     seeds(1,k) = seeds(1,k) + 1;
%                     break;
%                 end
%             end
% 
%         end
%     end
%     toc;
%     
%     [height,width]=size(saliencyMap);
%     for k = 1:nFrames
%         for i = 2:seeds(1,k)-1
%             currentIndex = seeds(i,k);
%             currentTrajectory = trajectories{currentIndex};
%             tr = currentTrajectory.trajectory;
%             trajectoryLength = size(tr,2);
%             trajectoryEnd = currentTrajectory.frameNum;
%             trajectoryStart = trajectoryEnd - trajectoryLength + 1;
%             hold on;
%             ylim([1 height]);
%             xlim([1 width]);
%             plot3(tr(2,:),(trajectoryStart-1):(trajectoryEnd-1),tr(1,:));
%         end
%     end
%     view(3);

%     % Group adjacent trajectories
%     groups = zeros(2000,GROUP_NUMBER);
%     groupIndex = 1;
%     groupItemIndex = ones(GROUP_NUMBER,1);
%     while 1
% 
%         if size(adjacencies,1) == 1; break; end
% 
%         % Find all adjacent trajectories with the current one
%         currentTrajectoryIndex = adjacencies(end,1);
%         indices = find((adjacencies == currentTrajectoryIndex));
%         [rowind,colind]=ind2sub(size(adjacencies),indices);
%         adjacentIndices = sub2ind(size(adjacencies),rowind,mod((colind),2)+1);
%         rowind = [rowind;rowind(1)];
%         colind = [colind;colind(1)];
%         adjacentIndices = [adjacentIndices;indices(1)];
% 
%         % Check if any of the adjacent trajectory is already member of a group
%         currentGroup = 0;
%         for i = 1:size(adjacentIndices)
%             currentGroup = find(groups == adjacencies(rowind(i),colind(i)));
%             if ~isempty(currentGroup)
%                 [~,currentGroup] = ind2sub(size(groups),currentGroup(1));
%                 break;
%             end
%         end
% 
%         % If non of the adjacent trajectories is a part of a group, create a
%         % new group for them
%         if isempty(currentGroup)
%             currentGroup = groupIndex;
%             groupIndex = groupIndex +1;
%         end
% 
%         % Assign trajectories to the defined group
%         currentGroupItemIndex = groupItemIndex(currentGroup);
%         groups(currentGroupItemIndex:currentGroupItemIndex +  ...
%             size(adjacentIndices,1)-1,currentGroup) = adjacencies(adjacentIndices);
%         groupItemIndex(currentGroup) = currentGroupItemIndex+size(adjacentIndices,1);
% 
%         % Remove the processed trajectories
%         adjacencies(rowind,:)=[];
% 
%     end
    groups = [];

end
