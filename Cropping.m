clc, clear,

RATIO = [16 9]; % aspect ratio of smartphones
CROP =10; %matrix having the crop value foe each frame
CROPPINGX = 20 * RATIO(1,1);
CROPPINGY = 20 * RATIO(1,2);
DOT = 0.1;
WINDOW_SIZE = 1;
FolderName = 'holywood2/fight/';
FileName = 'actioncliptrain00172.avi';

% UtilFunctions = Functions;
% UtilFunctions.ReadData(FolderName,FileName);
% UtilFunctions.CalculateMapping();
% load('Data.mat');
% 
% % Create movie structures for original, eye track info and cropped
% mov( 1 : floor( nFrames ))        = struct('cdata',zeros(vidHeight   ,vidWidth   , 3,'uint8'),'colormap',[]);
% movCenter( 1 : floor( nFrames ))  = struct('cdata',zeros(vidHeight   ,vidWidth   , 3,'uint8'),'colormap',[]);
% movCropped( 1 : floor( nFrames )) = struct('cdata',zeros(CROPPINGY*2 ,CROPPINGX*2, 3,'uint8'),'colormap',[]);
% 
% % Read movie in mov
% for k = 1 : nFrames
%     mov(k).cdata = read(video, k);    
% end
% 
% tic;
% [avgSaliency,distances] = UtilFunctions.CalculateMeanSaliency( nFrames , frames , eyes);
% % shotBoundaries = UtilFunctions.DetectShotBoundaries( mov );
% % [avgFlowOptical] = UtilFunctions.CreateFlow('optical',shotBoundaries, avgSaliency , mov , WINDOW_SIZE); % 'saliency' or 'optical'
% % [avgFlowSaliency] = UtilFunctions.CreateFlow('saliency',shotBoundaries, avgSaliency , mov , WINDOW_SIZE); % 'saliency' or 'optical'
% toc;

% save(strcat(FolderName , FileName,'.mat'));

load(strcat( FolderName, FileName,'.mat'));
% CROP =10; %matrix having the crop value foe each frame
% CROPPINGX = 20 * RATIO(1,1);
% CROPPINGY = 20 * RATIO(1,2);

% CROPARRAY = UtilFunctions.EstimateWindowSize(avgFlowOptical, eyes , shotBoundaries , CROP , RATIO , [vidWidth vidHeight]);
% [cropX cropY] = UtilFunctions.CreateWindow( CROPARRAY .* RATIO(1,1) , CROPARRAY .* RATIO(1,2) , avgFlowOptical(1:nFrames,:) , vidWidth , vidHeight );
% 
% 
% % Create cropped movie with cropX cropY values
% for k = 1 : nFrames
%     crop = mov(k).cdata( cropY(k,1):cropY(k,2) , cropX(k,1):cropX(k,2) , : );
%     movCropped(k).cdata = imresize(crop , [CROPPINGY*2 ,CROPPINGX*2]);
% end

% movie2avi(movCropped, strcat(FolderName , FileName,'_cropped.avi') , 'fps',  vidFPS);

% Create a video that shows the center of saliency in each frame
% DOT = repmat(5,nFrames,1);
% [ cropX cropY ] = UtilFunctions.CreateWindow( DOT , DOT , avgSaliency(2:nFrames+1,:) , vidWidth , vidHeight );
% for k = 1 : nFrames
%     movCenter(k).cdata = read(video, k);
%     movCenter(k).cdata( cropY(k,1):cropY(k,2) , cropX(k,1):cropX(k,2) , : ) = 254*ones() ;   
% end
% for k = 1 : nFrames 
%     imwrite(movCenter(k).cdata, sprintf('holywood2/fight/frames_center/video%05d.jpg', k)); 
% end 
% 
% DOT = repmat(0.5,nFrames,1);
% % Put a white dot for each center of attention in each frame
% for k = 1 : nFrames
%     
%     movCenter(k).cdata = read(video, k);
%     indices = find( frames( : , 1 ) == k );
%     nrIndices = length(indices);
%     
%     Fx = eyes( indices( 1:nrIndices ) , 2 );
%     Fy = eyes( indices( 1:nrIndices ) , 3 );
%     
%     [ cropX cropY ] = UtilFunctions.CreateWindow( DOT , DOT , [Fx Fy] , vidWidth , vidHeight );
%     for t = 1:nrIndices(1)
%           movCenter(k).cdata( cropY(t,1):cropY(t,2) , cropX(t,1):cropX(t,2) , : ) = 254*ones() ;   
%     end
%     [ cropX cropY ] = UtilFunctions.CreateWindow( DOT.*7 , DOT.*7, avgFlowOptical(k,:) , vidWidth , vidHeight );
%     test = cat(3,254*ones(cropX(1,2)-cropX(1,1)+1 , cropY(1,2)-cropY(1,1)+1) , zeros(cropX(1,2)-cropX(1,1)+1 , cropY(1,2)-cropY(1,1)+1), zeros(cropX(1,2)-cropX(1,1)+1 , cropY(1,2)-cropY(1,1)+1));
%     movCenter(k).cdata( cropY(1,1):cropY(1,2) , cropX(1,1):cropX(1,2) , : ) =test; 
% end

% % % Save created avi files
% movie2avi(movCropped, 'holywood2\cropped.avi' , 'fps',  vidFPS);
% movie2avi(movCenter, 'holywood2\others\center.avi' , 'fps',  vidFPS);
% movie2avi(mov, 'holywood2\original.avi' , 'fps',  vidFPS);
% for k = 1 : nFrames 
%     imwrite(mov(k).cdata, sprintf('holywood2/fight/frames_original_2/video%05d.jpg', k)); 
%     imwrite(movCropped(k).cdata, sprintf('holywood2/fight/frames_2/video%05d.jpg', k)); 
% end 

% Plotting the trajectory of center of attention
% UtilFunctions.PlotTrajectory2D( avgSaliency , 15);
% UtilFunctions.PlotTrajectory3D(avgSaliency);


% % Seam Carving and Linear Scaling for the Report
originalFrame = imread('C:\Users\ckocb_000\Desktop\courses\albert\project\Project\holywood2\fight\frames_original_2\video00035.jpg');
carved = im2uint8(seamcarving(originalFrame,150));
linScaled = imresize(originalFrame , [CROPPINGY CROPPINGX]);
imwrite(carved,'C:\Users\ckocb_000\Desktop\courses\albert\project\Project\holywood2\fight\frames_original_2\video00172_carved.jpg','jpg');
imwrite(linScaled,'C:\Users\ckocb_000\Desktop\courses\albert\project\Project\holywood2\fight\frames_original_2\video00172_scaled.jpg','jpg');




