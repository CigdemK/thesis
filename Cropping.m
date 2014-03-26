%% INITIALIZATION
clc, clear,

FolderName = '../holywood2/getoutofcar/';
FileName = 'actioncliptrain00121.avi';

RATIO = [16 9]; % aspect ratio of smartphones
CROP =10; %the crop value will be multiplied with aspect ration to get crop window

CROPPINGX = CROP * RATIO(1,1);
CROPPINGY = CROP * RATIO(1,2);

UtilFunctions = Functions;

%% CREATE OPTICAL FLOW OF THE MOVIE

UtilFunctions.ReadData(FolderName,FileName);
load('Data.mat');
mov = UtilFunctions.NewMovie(nFrames , vidHeight   ,vidWidth);
mov = UtilFunctions.ReadMovie(mov , video , nFrames );

% set optical flow parameters (see Coarse2FineTwoFrames.m for the definition of the parameters)
alpha = 0.012;
ratio = 0.75;
minWidth = 20;
nOuterFPIterations = 7;
nInnerFPIterations = 1;
nSORIterations = 30;

para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];

tic;
% Create cropped movie with cropX cropY values
movOpticalFlow = UtilFunctions.NewMovie(nFrames , vidHeight   ,vidWidth);
for k = 1 : nFrames-1
    
    im1 = mov(k).cdata;
    im2 = mov(k+1).cdata;
    [vx,vy,warpI2] = Coarse2FineTwoFrames(im1,im2,para);
    
    clear flow;
    flow(:,:,1) = vx;
    flow(:,:,2) = vy;
    imflow = flowToColor(flow);
    movOpticalFlow(k).cdata = imflow;
    
end
movOpticalFlow(k).cdata = imflow;
toc;

% Save the output
movie2avi(movOpticalFlow, strcat(FolderName , FileName,'_opticalFlow.avi') , 'fps',  vidFPS);

%% CREATING THE CROPPED MOVIE WITH A SALIENCY ALGORITHM
clear;clc;
FolderName = '../CAMO_Videos/Tilt/';
FileName = 'Tilt0015.avi';
UtilFunctions = Functions;
load(strcat(FolderName,FileName,'.mat'));
UtilFunctions.ReadData(FolderName,FileName);
load('Data.mat');
mov = UtilFunctions.NewMovie(nFrames , vidHeight   ,vidWidth);
mov = UtilFunctions.ReadMovie(mov , video , nFrames );

tic;
saliencyPoints = [];
allPoints = [];
for k = 1 : nFrames
    frame_map = gbvs(mov(k).cdata); 
    [r c] = find(frame_map.master_map_resized>0.2);
    frame_map = saliency(mov(k).cdata);
    frame_map = grayFrames(:,:,k);
    [r c] = find(frame_map>0.5);
    allPoints = cat(3,allPoints,frame_map);
    sizer = size(r);
    saliencyPoints = [saliencyPoints ; [repmat(k,sizer,1) r c]];
end
toc;

[ avgSaliency , distances ] = UtilFunctions.CalculateMeanSaliency( nFrames , saliencyPoints );
shotBoundaries = UtilFunctions.DetectShotBoundaries( mov );
[avgFlowOptical] = UtilFunctions.CreateFlow( 'optical' , shotBoundaries, avgSaliency , mov , 1,energyMap ); % 'saliency' or 'optical'
toc;
save(strcat(FolderName , FileName,'_2.mat'));
load(strcat(FolderName , FileName,'.mat'));UtilFunctions = Functions;

CROPARRAY = UtilFunctions.EstimateWindowSize(avgFlowOptical, saliencyPoints , shotBoundaries , CROP , RATIO , [vidWidth vidHeight]);
[cropX cropY] = UtilFunctions.CreateWindow( CROPARRAY .* RATIO(1,1) , CROPARRAY .* RATIO(1,2) , avgFlowOptical(1:nFrames,:) , vidWidth , vidHeight );
toc;
% Create cropped movie with cropX cropY values
movCropped = UtilFunctions.NewMovie(nFrames , vidHeight   ,vidWidth);
for k = 1 : nFrames
    crop = mov(k).cdata( cropY(k,1):cropY(k,2) , cropX(k,1):cropX(k,2) , : );
    movCropped(k).cdata = imresize(crop , [vidHeight   ,vidWidth]);
end

% Save the output
movie2avi(movCropped, strcat(FolderName , FileName,'_cropped1.avi') , 'fps',  vidFPS);
toc;

save(strcat(FolderName,FileName,'.mat'));
%% CREATING THE CROPPED MOVIE WITH EYETRACKING INFO

UtilFunctions.ReadData(FolderName,FileName);
UtilFunctions.ReadEyeTrackingData(FolderName,FileName);
UtilFunctions.CalculateMapping();
load('Data.mat');

% Create movie structure & find optical flow per shot
mov = UtilFunctions.ReadMovie( video , nFrames , vidHeight , vidWidth);
[avgSaliency,distances] = UtilFunctions.CalculateMeanSaliency( nFrames , eyes);
shotBoundaries = UtilFunctions.DetectShotBoundaries( mov );
[avgFlowOptical] = UtilFunctions.CreateFlow( 'optical' , shotBoundaries, avgSaliency , mov , 1 ); % 'saliency' or 'optical'
[avgFlowSaliency] = UtilFunctions.CreateFlow( 'saliency' , shotBoundaries, avgSaliency , mov , 1 ); % 'saliency' or 'optical'

save(strcat(FolderName , FileName,'.mat'));
% load(strcat( FolderName, FileName,'.mat'));

CROPARRAY = UtilFunctions.EstimateWindowSize(avgFlowOptical, eyes , shotBoundaries , CROP , RATIO , [vidWidth vidHeight]);
[cropX cropY] = UtilFunctions.CreateWindow( CROPARRAY .* RATIO(1,1) , CROPARRAY .* RATIO(1,2) , avgFlowOptical(1:nFrames,:) , vidWidth , vidHeight );

% Create cropped movie with cropX cropY values
for k = 1 : nFrames
    crop = mov(k).cdata( cropY(k,1):cropY(k,2) , cropX(k,1):cropX(k,2) , : );
    movCropped(k).cdata = imresize(crop , [CROPPINGY*2 ,CROPPINGX*2]);
end

% Save the output
movie2avi(movCropped, strcat(FolderName , FileName,'_cropped.avi') , 'fps',  vidFPS);


%% CREATING MOVIE THAT SHOWS CENTER OF SALIENCY PER FRAME

DOT = repmat(5,nFrames,1);
movCenter = UtilFunctions.ReadMovie( video , nFrames , vidHeight , vidWidth);

[ cropX cropY ] = UtilFunctions.CreateWindow( DOT , DOT , avgSaliency(2:nFrames+1,:) , vidWidth , vidHeight );

for k = 1 : nFrames
    movCenter(k).cdata( cropY(k,1):cropY(k,2) , cropX(k,1):cropX(k,2) , : ) = 254*ones() ;   
end

% Save the output
movie2avi(movCenter, strcat(FolderName , FileName,'_center.avi') , 'fps',  vidFPS);


%% CREATING MOVIE THAT SHOWS ALL SALIENCY POINTS PER FRAME

DOT = repmat(0.5,nFrames,1);
movAllPoints = UtilFunctions.ReadMovie( video , nFrames , vidHeight , vidWidth);

for k = 1 : nFrames
    
    indices = find( frames( : , 1 ) == k );
    nrIndices = length(indices);
    
    Fx = eyes( indices( 1:nrIndices ) , 2 );
    Fy = eyes( indices( 1:nrIndices ) , 3 );
    
    [ cropX cropY ] = UtilFunctions.CreateWindow( DOT , DOT , [Fx Fy] , vidWidth , vidHeight );
    
    for t = 1:nrIndices(1)
          movAllPoints(k).cdata( cropY(t,1):cropY(t,2) , cropX(t,1):cropX(t,2) , : ) = 254*ones() ;   
    end
    
%     Uncomment if you want to add a big red dot for center of saliency
%     [ cropX cropY ] = UtilFunctions.CreateWindow( DOT.*7 , DOT.*7, avgFlowOptical(k,:) , vidWidth , vidHeight );
%     test = cat(3,254*ones(cropX(1,2)-cropX(1,1)+1 , cropY(1,2)-cropY(1,1)+1) , zeros(cropX(1,2)-cropX(1,1)+1 , cropY(1,2)-cropY(1,1)+1), zeros(cropX(1,2)-cropX(1,1)+1 , cropY(1,2)-cropY(1,1)+1));
%     movAllPoints(k).cdata( cropY(1,1):cropY(1,2) , cropX(1,1):cropX(1,2) , : ) = test; 
    
end

% Save the output
movie2avi(movAllPoints, strcat(FolderName , FileName,'_allpoints.avi') , 'fps',  vidFPS);


%% SAVING FRAMES SEPERATELY

for k = 1 : nFrames 
    imwrite(movCenter(k).cdata, sprintf('%s%s/frames_center/frames%05d.jpg', FolderName , FileName, k)); 
end 

%% FOR REPORT

% Plotting the trajectory of center of attention
UtilFunctions.PlotTrajectory2D( avgSaliency , 15);
UtilFunctions.PlotTrajectory3D( avgSaliency );

% Seam Carving and Linear Scaling for the Report
originalFrame = imread('C:\Users\ckocb_000\Desktop\courses\albert\project\Project\holywood2\fight\frames_original_2\video00035.jpg');
carved = im2uint8(seamcarving(originalFrame,150));
linScaled = imresize(originalFrame , [CROPPINGY CROPPINGX]);
imwrite(carved,'C:\Users\ckocb_000\Desktop\courses\albert\project\Project\holywood2\fight\frames_original_2\video00172_carved.jpg','jpg');
imwrite(linScaled,'C:\Users\ckocb_000\Desktop\courses\albert\project\Project\holywood2\fight\frames_original_2\video00172_scaled.jpg','jpg');


