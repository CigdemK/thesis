function F = Cropping
    F.CropWithHomography = @CropWithHomography;
    F.ShowSaliencyPoints = @ShowSaliencyPoints;
    F.CropWithVideoSaliency = @CropWithVideoSaliency;
end

function ShowSaliencyPoints(FolderName , FileName , mode) %mode = 'all' or 'mean'

    Reader = ReadFunctions;
    Cropper = CropFunctions;
    Video = VideoFunctions;
    
    [video,nFrames,vidHeight,vidWidth,vidFPS] = Reader.ReadData(FolderName,FileName);
    
    movAllPoints = Reader.NewMovie(nFrames , vidHeight   ,vidWidth);
    movAllPoints = Reader.ReadMovie(movAllPoints , video );
    
    Video = VideoFunctions;
    saliencyPoints = Video.CalculateVideoSaliency(movAllPoints);
%     load(strcat(FolderName , FileName,'.mat') , 'saliencyPoints' );

    DOT = repmat(0.5,nFrames,1);
    if strcmp( mode , 'all' )
        
        for k = 1 : nFrames

            indices = find( saliencyPoints( : , 1 ) == k );
            nrIndices = length(indices);

            Fx = saliencyPoints( indices( 1:nrIndices ) , 2 );
            Fy = saliencyPoints( indices( 1:nrIndices ) , 3 );

            for i = 1:nrIndices
                [ cropX cropY ] = Cropper.CreateWindow( DOT , DOT , [Fx(i) Fy(i)] , vidWidth , vidHeight );
                movAllPoints(k).cdata( cropY(1):cropY(2) , cropX(1):cropX(2) , : ) = 254*ones() ;   
            end
            toc;
        end
        
    elseif strcmp( mode , 'mean' )
        
        for k = 1 : nFrames
            
            avgSaliency = Video.CalculateMeanSaliency( nFrames , saliencyPoints );
            [ cropX cropY ] = Cropper.CreateWindow( DOT.*7 , DOT.*7, avgSaliency , vidWidth , vidHeight );
            redDot = cat( 3 , 254*ones(cropX(1,2)-cropX(1,1)+1 , cropY(1,2)-cropY(1,1)+1) , ...
                            zeros(cropX(1,2)-cropX(1,1)+1 , cropY(1,2)-cropY(1,1)+1), ...
                            zeros(cropX(1,2)-cropX(1,1)+1 , cropY(1,2)-cropY(1,1)+1));
            movAllPoints(k).cdata( cropY(1,1):cropY(1,2) , cropX(1,1):cropX(1,2) , : ) = redDot; 
        
        end

    end

    % Save the output
    movie2avi(movAllPoints, strcat(FolderName , FileName,'_allpoints_video.avi') , 'fps',  vidFPS);

end

function CropWithHomography(FolderName, FileName , CROP , RATIO , mode) %mode = 'optical', 'saliency' or 'homography'

    Reader = ReadFunctions;
    [video,nFrames,vidHeight,vidWidth,vidFPS] = Reader.ReadData(FolderName,FileName);
    fprintf('Number of frames = %d\n',nFrames);

    % Create movie structure & find optical flow per shot
    mov = Reader.NewMovie(nFrames , vidHeight   ,vidWidth);
    mov = Reader.ReadMovie(mov , video );

    Video = VideoFunctions;
    saliencyPoints = Video.CalculateStaticSaliency(mov);
    avgSaliency = Video.CalculateMeanSaliency( nFrames , saliencyPoints );
    shotBoundaries = Video.DetectShotBoundaries( mov );
    [avgFlowOptical_h] = Video.CreateFlow( mode , shotBoundaries , avgSaliency , mov , 4 ); % 'saliency' or 'optical' or 'homography'

%     save(strcat(FolderName , FileName,'.mat'));
%     load(strcat(FolderName , FileName,'.mat'));

    Cropper = CropFunctions;
    cropArray = Cropper.EstimateWindowSize(avgFlowOptical_h, saliencyPoints , shotBoundaries , CROP , RATIO , [vidWidth vidHeight]);
    [cropX cropY] = Cropper.CreateWindow( cropArray .* RATIO(1,1) , cropArray .* RATIO(1,2) , avgFlowOptical_h(1:nFrames,:) , vidWidth , vidHeight );
    
    % Create cropped movie with cropX cropY values
    movCropped = Reader.NewMovie(nFrames , vidHeight   ,vidWidth);
    for k = 1 : nFrames
        crop = mov(k).cdata( cropY(k,1):cropY(k,2) , cropX(k,1):cropX(k,2) , : );
        movCropped(k).cdata = imresize(crop , [CROP.*RATIO(1,2)*2 ,CROP.*RATIO(1,1)*2]);
    end

    % Save the output
    movie2avi(movCropped, strcat(FolderName , FileName,'_cropped_optical.avi') , 'fps',  vidFPS);

end

function CropWithVideoSaliency(FolderName,FileName,CROP,RATIO,mode) %mode = 'optical', 'saliency' or 'homography'

    Reader = ReadFunctions;
    [video,nFrames,vidHeight,vidWidth,vidFPS] = Reader.ReadData(FolderName,FileName);
    fprintf('Number of frames = %d\n',nFrames);

    % Create movie structure & find optical flow per shot
    mov = Reader.NewMovie(nFrames , vidHeight   ,vidWidth);
    mov = Reader.ReadMovie(mov , video );

    Video = VideoFunctions;    
    saliencyPoints = Video.CalculateVideoSaliency(mov);
    avgSaliency = Video.CalculateMeanSaliency( nFrames , saliencyPoints );   
    shotBoundaries = Video.DetectShotBoundaries( mov );
    [avgFlowOptical_h] = Video.CreateFlow( mode , shotBoundaries , avgSaliency , mov , 4 ); % 'saliency' or 'optical' or 'homography'

%     save(strcat(FolderName , FileName,'.mat'));
%     load(strcat(FolderName , FileName,'.mat'));

    Cropper = CropFunctions;
    cropArray = Cropper.EstimateWindowSize(avgFlowOptical_h, saliencyPoints , shotBoundaries , CROP , RATIO , [vidWidth vidHeight]);
    [cropX cropY] = Cropper.CreateWindow( cropArray .* RATIO(1,1) , cropArray .* RATIO(1,2) , avgFlowOptical_h(1:nFrames,:) , vidWidth , vidHeight );
    
    % Create cropped movie with cropX cropY values
    movCropped = Reader.NewMovie(nFrames , vidHeight   ,vidWidth);
    for k = 1 : nFrames
        crop = mov(k).cdata( cropY(k,1):cropY(k,2) , cropX(k,1):cropX(k,2) , : );
        movCropped(k).cdata = imresize(crop , [CROP.*RATIO(1,2)*2 ,CROP.*RATIO(1,1)*2]);
    end

    % Save the output
    movie2avi(movCropped, strcat(FolderName , FileName,'_cropped_optical.avi') , 'fps',  vidFPS);

end
% 
% %% CREATE OPTICAL FLOW OF THE MOVIE
% 
% UtilFunctions.ReadData(FolderName,FileName);
% load('Data.mat');
% mov = UtilFunctions.NewMovie(nFrames , vidHeight   ,vidWidth);
% mov = UtilFunctions.ReadMovie(mov , video , nFrames );
% 
% % set optical flow parameters (see Coarse2FineTwoFrames.m for the definition of the parameters)
% alpha = 0.012;
% ratio = 0.75;
% minWidth = 20;
% nOuterFPIterations = 7;
% nInnerFPIterations = 1;
% nSORIterations = 30;
% 
% para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];
% 
% tic;
% % Create cropped movie with cropX cropY values
% movOpticalFlow = UtilFunctions.NewMovie(nFrames , vidHeight   ,vidWidth);
% for k = 1 : nFrames-1
%     
%     im1 = mov(k).cdata;
%     im2 = mov(k+1).cdata;
%     [vx,vy,warpI2] = Coarse2FineTwoFrames(im1,im2,para);
%     
%     clear flow;
%     flow(:,:,1) = vx;
%     flow(:,:,2) = vy;
%     imflow = flowToColor(flow);
%     movOpticalFlow(k).cdata = imflow;
%     
% end
% movOpticalFlow(k).cdata = imflow;
% toc;
% 
% % Save the output
% movie2avi(movOpticalFlow, strcat(FolderName , FileName,'_opticalFlow.avi') , 'fps',  vidFPS);
% 
% %% CREATING THE CROPPED MOVIE WITH A SALIENCY ALGORITHM
% clear;clc;
% FolderName = '../CAMO_Videos/Tilt/';
% FileName = 'Tilt0015.avi';
% UtilFunctions = Functions;
% load(strcat(FolderName,FileName,'.mat'));
% UtilFunctions.ReadData(FolderName,FileName);
% load('Data.mat');
% mov = UtilFunctions.NewMovie(nFrames , vidHeight   ,vidWidth);
% mov = UtilFunctions.ReadMovie(mov , video , nFrames );
% 
% tic;
% saliencyPoints = [];
% allPoints = [];
% for k = 1 : nFrames
%     frame_map = gbvs(mov(k).cdata); 
%     [r c] = find(frame_map.master_map_resized>0.2);
%     frame_map = saliency(mov(k).cdata);
%     frame_map = grayFrames(:,:,k);
%     [r c] = find(frame_map>0.5);
%     allPoints = cat(3,allPoints,frame_map);
%     sizer = size(r);
%     saliencyPoints = [saliencyPoints ; [repmat(k,sizer,1) r c]];
% end
% toc;
% 
% [ avgSaliency , distances ] = UtilFunctions.CalculateMeanSaliency( nFrames , saliencyPoints );
% shotBoundaries = UtilFunctions.DetectShotBoundaries( mov );
% [avgFlowOptical] = UtilFunctions.CreateFlow( 'optical' , shotBoundaries, avgSaliency , mov , 1,energyMap ); % 'saliency' or 'optical'
% toc;
% save(strcat(FolderName , FileName,'_2.mat'));
% load(strcat(FolderName , FileName,'.mat'));UtilFunctions = Functions;
% 
% CROPARRAY = UtilFunctions.EstimateWindowSize(avgFlowOptical, saliencyPoints , shotBoundaries , CROP , RATIO , [vidWidth vidHeight]);
% [cropX cropY] = UtilFunctions.CreateWindow( CROPARRAY .* RATIO(1,1) , CROPARRAY .* RATIO(1,2) , avgFlowOptical(1:nFrames,:) , vidWidth , vidHeight );
% toc;
% % Create cropped movie with cropX cropY values
% movCropped = UtilFunctions.NewMovie(nFrames , vidHeight   ,vidWidth);
% for k = 1 : nFrames
%     crop = mov(k).cdata( cropY(k,1):cropY(k,2) , cropX(k,1):cropX(k,2) , : );
%     movCropped(k).cdata = imresize(crop , [vidHeight   ,vidWidth]);
% end
% 
% % Save the output
% movie2avi(movCropped, strcat(FolderName , FileName,'_cropped1.avi') , 'fps',  vidFPS);
% toc;
% 
% save(strcat(FolderName,FileName,'.mat'));
% %% CREATING THE CROPPED MOVIE WITH EYETRACKING INFO
% 
% UtilFunctions.ReadData(FolderName,FileName);
% UtilFunctions.ReadEyeTrackingData(FolderName,FileName);
% UtilFunctions.CalculateMapping();
% load('Data.mat');
% 
% % Create movie structure & find optical flow per shot
% mov = UtilFunctions.ReadMovie( video , nFrames , vidHeight , vidWidth);
% [avgSaliency,distances] = UtilFunctions.CalculateMeanSaliency( nFrames , eyes);
% shotBoundaries = UtilFunctions.DetectShotBoundaries( mov );
% [avgFlowOptical] = UtilFunctions.CreateFlow( 'optical' , shotBoundaries, avgSaliency , mov , 1 ); % 'saliency' or 'optical'
% [avgFlowSaliency] = UtilFunctions.CreateFlow( 'saliency' , shotBoundaries, avgSaliency , mov , 1 ); % 'saliency' or 'optical'
% 
% save(strcat(FolderName , FileName,'.mat'));
% % load(strcat( FolderName, FileName,'.mat'));
% 
% CROPARRAY = UtilFunctions.EstimateWindowSize(avgFlowOptical, eyes , shotBoundaries , CROP , RATIO , [vidWidth vidHeight]);
% [cropX cropY] = UtilFunctions.CreateWindow( CROPARRAY .* RATIO(1,1) , CROPARRAY .* RATIO(1,2) , avgFlowOptical(1:nFrames,:) , vidWidth , vidHeight );
% 
% % Create cropped movie with cropX cropY values
% for k = 1 : nFrames
%     crop = mov(k).cdata( cropY(k,1):cropY(k,2) , cropX(k,1):cropX(k,2) , : );
%     movCropped(k).cdata = imresize(crop , [CROPPINGY*2 ,CROPPINGX*2]);
% end
% 
% % Save the output
% movie2avi(movCropped, strcat(FolderName , FileName,'_cropped.avi') , 'fps',  vidFPS);
% 
% 
% %% CREATING MOVIE THAT SHOWS CENTER OF SALIENCY PER FRAME
% 
% DOT = repmat(5,nFrames,1);
% movCenter = UtilFunctions.ReadMovie( video , nFrames , vidHeight , vidWidth);
% 
% [ cropX cropY ] = UtilFunctions.CreateWindow( DOT , DOT , avgSaliency(2:nFrames+1,:) , vidWidth , vidHeight );
% 
% for k = 1 : nFrames
%     movCenter(k).cdata( cropY(k,1):cropY(k,2) , cropX(k,1):cropX(k,2) , : ) = 254*ones() ;   
% end
% 
% % Save the output
% movie2avi(movCenter, strcat(FolderName , FileName,'_center.avi') , 'fps',  vidFPS);
% 
% 

% %% SAVING FRAMES SEPERATELY
% 
% for k = 1 : nFrames 
%     imwrite(movCenter(k).cdata, sprintf('%s%s/frames_center/frames%05d.jpg', FolderName , FileName, k)); 
% end 
% 
% %% FOR REPORT
% 
% % Plotting the trajectory of center of attention
% UtilFunctions.PlotTrajectory2D( avgSaliency , 15);
% UtilFunctions.PlotTrajectory3D( avgSaliency );
% 
% % Seam Carving and Linear Scaling for the Report
% originalFrame = imread('C:\Users\ckocb_000\Desktop\courses\albert\project\Project\holywood2\fight\frames_original_2\video00035.jpg');
% carved = im2uint8(seamcarving(originalFrame,150));
% linScaled = imresize(originalFrame , [CROPPINGY CROPPINGX]);
% imwrite(carved,'C:\Users\ckocb_000\Desktop\courses\albert\project\Project\holywood2\fight\frames_original_2\video00172_carved.jpg','jpg');
% imwrite(linScaled,'C:\Users\ckocb_000\Desktop\courses\albert\project\Project\holywood2\fight\frames_original_2\video00172_scaled.jpg','jpg');
% 


