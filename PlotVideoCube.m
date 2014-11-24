clc,clear;
folderName = 'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\';
fileName = 'actioncliptest00001.avi'; %_ImprovedTrajectoriesSaliencyMap.mat
shotStart = 48;
shotEnd = 95;

video=VideoReader([folderName fileName]);
frames = read(video);
frames = frames(:,:,:,shotStart:shotEnd);
[height,width,~,nFrames] = size(frames);

tic;
regexres = regexp(fileName,'.avi','split');
load([folderName regexres{1} '\ImprovedTrajectoryOriginalTJS.mat'])
toc;

% Try grouping neighbour trajectories
% accordingToStart = cell(nFrames,1);
% for i = 1:nFrames; accordingToStart{i} = []; end
forKMeans = zeros(size(trajectories,1),5); %[meanX,meanY,firstD,secondD]
correspondance = zeros(size(trajectories,1),2);
kMeansIndex = 1;
for k = 1:size(trajectories,1)
    
    currentTrajectory = trajectories{k};
    
    trajectoryLength = size(currentTrajectory.trajectory,2);
    trajectoryEnd = currentTrajectory.frameNum;
    trajectoryStart = trajectoryEnd - trajectoryLength + 1;

    if trajectoryStart > shotStart && trajectoryEnd < shotEnd   
        tr = currentTrajectory.trajectory;
        forKMeans(k,:) = dtw_c(tr(1,:),tr(2,:));
%         X = fft(tr(1,:));
%         Y = fft(tr(2,:));
     
%         currentTrajectory.distance = sqrt(tr(1,:).^2 + tr(2,:).^2);
%         meanDistance = mean(currentTrajectory.distance);
%         p = polyfit(1:trajectoryLength,currentTrajectory.distance,3);
%         pd = polyder(p);
%         pdd = polyder(pd);
            
%         forKMeans(kMeansIndex,:) = [real(X(1)),real(Y(1))];  
%         forKMeans(kMeansIndex,:) = [pd,pdd]; 
%         kMeansIndex = kMeansIndex + 1;
%         correspondance(kMeansIndex,:) = [kMeansIndex,k];
%     else
%         forKMeans(k,:) = [0,0,0,0,0,0,0,0];
    end
    
%     if trajectoryStart > shotStart && trajectoryEnd < shotEnd
%         accordingToStart{trajectoryStart - shotStart}(end+1) = k;
%     end

end
forKMeans(~any(forKMeans,2),:) = [];
correspondance(~any(correspondance,2),:) = [];
toc;


GROUP_NUMBER = 15;
params = kmeans(forKMeans,GROUP_NUMBER);
groups = params(end).classes;
toc;

% Plot trajectory groups in different colors
colorArr = {'b','g','k','r','m','c','y'};   
for k = 1:50:size(trajectories,1)

    currentTrajectory = trajectories{k};
    tr = currentTrajectory.trajectory;
    trajectoryLength = size(currentTrajectory.trajectory,2);
    trajectoryEnd = currentTrajectory.frameNum;
    trajectoryStart = trajectoryEnd - trajectoryLength + 1;
    
    if trajectoryStart < shotStart || trajectoryEnd > shotEnd; continue; end
    
    trGroup = groups(k);
    hold on;
%     mod(trGroup,7)+1
    plot3(tr(2,:),(trajectoryStart-1):(trajectoryEnd-1),tr(1,:),[colorArr{mod(trGroup,7)+1}]);
 
end
axis tight;
view(3);
toc;

% ADJACENCY_DISTANCE = 20;
% adjacencies = zeros(100000,2);
% adjacencies(1,1) = 2; % keep the last row
% for k = 1:size(accordingToStart,1)
%     currentStart = accordingToStart{k};
%     for i = 1:97:size(currentStart,2)
%         
%         traji = trajectories{currentStart(i)}.trajectory; 
%         
%         % Chech trajectories starting from the same frame
%         for j = i+1:size(currentStart,2)
%             trajj = trajectories{currentStart(j)}.trajectory;
%             distance = abs(mean(traji(2,:) - trajj(2,:)));
%             if distance < ADJACENCY_DISTANCE
%                 adjacencies(adjacencies(1,1),:) = [currentStart(i),currentStart(j)];
%                 adjacencies(1,1) = adjacencies(1,1) + 1;
%             end
%         end
% 
%         % Chech trajectories starting from the -1 frame
%         if k > 1
%             possibleStart = accordingToStart{k-1};
%             for t = 1:size(possibleStart,2)
%                 trajj = trajectories{possibleStart(t)}.trajectory;
%                 distance = abs(mean(traji(2,1:end-1) - trajj(2,2:end)));
%                 if distance < ADJACENCY_DISTANCE
%                     adjacencies(adjacencies(1,1),:) = [currentStart(i),possibleStart(t)];
%                     adjacencies(1,1) = adjacencies(1,1) + 1;
%                 end
%             end
%         end
%     end
% %     toc;
% end
% toc;
%   
% adjacencies(~any(adjacencies,2),:) = []; % remove zero rows
% 
% % Group adjacent trajectories
% GROUP_NUMBER = 1000;
% groups = zeros(2000,GROUP_NUMBER);
% groupIndex = 1;
% groupItemIndex = ones(GROUP_NUMBER,1);
% while 1
% 
%     if size(adjacencies,1) == 1; break; end
% 
%     % Find all adjacent trajectories with the current one
%     currentTrajectoryIndex = adjacencies(end,1);
%     indices = find((adjacencies == currentTrajectoryIndex));
%     [rowind,colind]=ind2sub(size(adjacencies),indices);
%     adjacentIndices = sub2ind(size(adjacencies),rowind,mod((colind),2)+1);
%     rowind = [rowind;rowind(1)];
%     colind = [colind;colind(1)];
%     adjacentIndices = [adjacentIndices;indices(1)];
% 
%     % Check if any of the adjacent trajectory is already member of a group
%     currentGroup = 0;
%     for i = 1:size(adjacentIndices)
%         currentGroup = find(groups == adjacencies(rowind(i),colind(i)));
%         if ~isempty(currentGroup)
%             [~,currentGroup] = ind2sub(size(groups),currentGroup(1));
%             break;
%         end
%     end
% 
%     % If non of the adjacent trajectories is a part of a group, create a
%     % new group for them
%     if isempty(currentGroup)
%         currentGroup = groupIndex;
%         groupIndex = groupIndex +1;
%     end
% 
%     % Assign trajectories to the defined group
%     currentGroupItemIndex = groupItemIndex(currentGroup);
%     groups(currentGroupItemIndex:currentGroupItemIndex +  ...
%         size(adjacentIndices,1)-1,currentGroup) = adjacencies(adjacentIndices);
%     groupItemIndex(currentGroup) = currentGroupItemIndex+size(adjacentIndices,1);
% 
%     % Remove the processed trajectories
%     adjacencies(rowind,:)=[];
% 
% end
% toc;
% 
% 
% % Remove all zero columns
% groups(:,~any(groups,1)) = []; 
% toc;
% % Plot trajectory groups in different colors
% colorArr = {'b','g','k','r','m','c','y'};   
% for k = 1:size(groups,2)
% 
%     currentIndices = groups(:,k);
%     currentIndices(~any(currentIndices,2),:) = [];  %remove zero rows
%     if size(currentIndices,1) < 1000; continue; end
%     
%     for i = 1:100: size(currentIndices)
%         currentTrajectory = trajectories{currentIndices(i)};
%         tr = currentTrajectory.trajectory;
%         trajectoryLength = size(tr,2);
%         trajectoryEnd = currentTrajectory.frameNum;
%         trajectoryStart = trajectoryEnd - trajectoryLength + 1;
%         hold on;
%         plot3(tr(2,:),(trajectoryStart-1):(trajectoryEnd-1),tr(1,:),[colorArr{mod(k,7)+1}]);
%     end
%     
% end
% % axis equal;
% view(3);
% 
% % Plot video cube
% figure;
% imgSaliency = ImageSaliency;
% saliencyMapFixations = imgSaliency.ActionsInTheEyeFixation([folderName fileName]);
% for i = shotStart:shotEnd
%     cdata = flipdim( saliencyMapFixations(:,:,i), 1 );
%     surface([0 width; 0 width], [(i-1)*10 (i-1)*10; (i-1)*10 (i-1)*10], [0 0; height height], ...
%         'FaceColor', 'texturemap', 'CData', cdata , 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
%     xlabel('Video Width','FontWeight','Bold');
%     ylabel('Time / Frames','FontWeight','Bold');
%     zlabel('Video Height','FontWeight','Bold');
%     axis equal;
%     view(3);
% end
% 
% Plot trajectories
% for k = 1:100:size(trajectories,1)
% 
%     currentTrajectory = trajectories{k};
%     tr = currentTrajectory.trajectory;
%    
%     trajectoryLength = size(tr,2);
% 
%     trajectoryEnd = currentTrajectory.frameNum;
%     trajectoryStart = trajectoryEnd - trajectoryLength + 1;
%     if (trajectoryEnd > shotEnd) || (trajectoryStart < shotStart); continue; end
%     
%     hold on;
% %     plot3(tr(2,:),(trajectoryStart-1)*10:10:(trajectoryEnd-1)*10,tr(1,:),'r')
%     plot(tr(2,:),(trajectoryStart-1)*10:10:(trajectoryEnd-1)*10,'x')
% end
% 
% % Plot trajectory saliency and fixation saliency
% % imgSaliency = ImageSaliency;
% % saliencyMap = imgSaliency.Cigdem([folderName fileName],trajectoriesByFrame);
% saliencyMapFixations = imgSaliency.ActionsInTheEyeFixation([folderName fileName]);
% for i = shotStart:5:shotEnd
%     figure;
%     subplot(2,1,1);
%     imshow(frames(:,:,:,i)); 
%     subplot(2,1,2);
%     imshow(saliencyMapFixations(:,:,i));
%     subplot(4,1,3);
%     imshow(saliencyMap(:,:,i));
%     subplot(4,1,4);
%     
% 
%     cascade='haarcascade_frontalface_alt2.xml'; % a little noisy a few misses
%     FaceData = FaceDetect(cascade,double(rgb2gray(frames(:,:,:,i)))); 
% 
%     if find(FaceData==0)
%         FaceData(find(FaceData==0))=1;
%     %     fprintf('changingFaceData from 0');
%     end
% 
%     if FaceData ~= -1
%         for j=1:size(FaceData, 1)
%             Face = FaceData(j, :);
%             Faces(Face(2):Face(2)+Face(4), Face(1):Face(1)+Face(3)) = 1;
%         end
%     end
%     Faces = imresize(Faces, [height,width]);
%     imshow(Faces);
% end
