function F = Ali59MFunctions
    F.ShowSaliencyPoints = @ShowSaliencyPoints;
%     F.CropWithHomography = @CropWithHomography;
    F.CropWithVideoSaliency = @CropWithVideoSaliency;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Public Functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function frames = ShowSaliencyPoints(moviePath , mode, saliency) 
%mode = 'all' or 'mean', saliency = 'gaze' or 'estimate'   

    if nargin < 3
        saliency = 'estimate';
    end
    if nargin < 2
        mode = 'all';
    end
    
    video = VideoReader( moviePath );
    frames = read(video);
    
    VideoSal = VideoSaliency;
    AITE = ActionsInTheEye;
    if strcmp( saliency , 'estimate' )
        [saliencyPoints,~] = VideoSal.Nguyen2013('',frames);
    elseif strcmp( saliency , 'gaze' )
        resultMap = AITE.ReadEyeTrackingData(moviePath);
        resultMap.vidDuration = video.Duration;
        resultMap.nFrames = video.NumberOfFrames;
        saliencyPoints = AITE.CalculateMapping(resultMap);
    end

    DOT = repmat(0.5,video.NumberOfFrames,1);
    Cropper = CropFunctions;
    if strcmp( mode , 'all' )
        
        for k = 1 : video.NumberOfFrames

            indices = find( saliencyPoints( : , 2 ) == k );
            nrIndices = length(indices);

            Fx = saliencyPoints( indices , 4 );
            Fy = saliencyPoints( indices , 3 );

            for i = 1:nrIndices
                [cropX, cropY] = Cropper.CreateCropWindow(DOT ,DOT ,[Fx(i) Fy(i)] ,video.Width ,video.Height);
                frames(cropY(1):cropY(2) ,cropX(1):cropX(2) ,:, k) = 254*ones() ;   
            end
      
        end
        
    elseif strcmp( mode , 'mean' )
        
        avgSaliency = CalculateMeanSaliency(video.NumberOfFrames ,saliencyPoints);
        for k = 1 : video.NumberOfFrames
            
            Fx = avgSaliency(k,2);
            Fy = avgSaliency(k,1);
            [cropX, cropY] = Cropper.CreateCropWindow(DOT.*7 ,DOT.*7,[Fx , Fy] , video.Width , video.Height );
            redDot = cat( 3 , 254 * ones(cropX(1,2)-cropX(1,1)+1 , cropY(1,2)-cropY(1,1)+1), ...
                                   zeros(cropX(1,2)-cropX(1,1)+1 , cropY(1,2)-cropY(1,1)+1), ...
                                   zeros(cropX(1,2)-cropX(1,1)+1 , cropY(1,2)-cropY(1,1)+1));
            frames(cropY(1):cropY(2), cropX(1):cropX(2), :, k) = redDot; 
        
        end

    end

end

function frames = CropWithVideoSaliency(moviePath,CROP)

    video = VideoReader( moviePath );
    frames = read(video);
frames = frames(:,:,:,1:3);
    
    Cropper = CropFunctions; 
    VideoSal = VideoSaliency;
    Util = UtilFunctions;
    
    frames = Cropper.RemoveBlackBars(frames);
    [vidHeight,vidWidth,~,nFrames] = size(frames);
    shotBoundaries = Util.DetectShotBoundaries(frames);
    
    [videoSaliencyMap , opticalFlowMap] = VideoSal.Nguyen2013('',frames);
    importantPts = ThresholdSaliency(videoSaliencyMap);
    avgSaliency = CalculateMeanSaliency(nFrames , importantPts );   
    [avgOpticalFlow] = CreateFlow(shotBoundaries, avgSaliency, frames, opticalFlowMap);

    cropArray = Cropper.EstimateCropWindowSize(avgOpticalFlow, importantPts, shotBoundaries, CROP);
    [cropX, cropY] = Cropper.CreateCropWindow( cropArray .* CROP(1,1) , cropArray .* CROP(1,2) , avgOpticalFlow(1:nFrames,:) , vidWidth , vidHeight );

    % Create cropped movie with cropX cropY values
    for k = 1 : nFrames
        crop = frames(cropY(k,1):cropY(k,2), cropX(k,1):cropX(k,2), :, k);
        frames(:,:,:,k) = imresize(crop , [CROP ,CROP]);
    end

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Private Functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [avg] = CalculateMeanSaliency( maxFrame  , eyes )
    avg = [];
    for k = 1 : maxFrame
        indices = find( eyes( : , 2 ) == k );
        avg = [ avg ; mean( eyes ( indices(:) , 3:4 ) ) ];
    end   
end

function [avgFlow] = CreateFlow(shotBoundaries, avgSaliency , frames , opticalFlowMap )

    if nargin < 6
        Motion = MotionFunctions;
        [opticalFlowX,opticalFlowY] = Motion.MyOpticalFlow(frames);
        opticalFlowMap(:,:,:,1) = opticalFlowX;
        opticalFlowMap(:,:,:,2) = opticalFlowY;
    end
    
    avgFlow = [];
    nrShots = size(shotBoundaries,2)-1;
    
    Motion = MotionFunctions;
    for k = 1:nrShots
        shotStart = shotBoundaries(k);
        shotEnd = shotBoundaries(k+1)-1;

        shotSaliency = avgSaliency( shotStart:shotEnd , :);

        [Vx,Vy] = Motion.CameraSpeed(opticalFlowMap);
        avgFlowTmp = CalculateFlowMean( shotSaliency(1,:) , Vx ,Vy);   

        avgFlow = [avgFlow;avgFlowTmp];
    end

end

function [avgFlow] = CalculateFlowMean(startingPoint , Vx ,Vy)

    shotLength = length(Vx);
    Y = startingPoint(1,1);
    X = startingPoint(1,2);
    avgFlow = [X Y];
    meanFlow = [Vx Vy];
    
    for k = 2 : shotLength %dismiss the first frame of the shot since it is set to startingPoint
        transformed = [meanFlow(k,:) 1]';
        tmp = transformed(1:2)' + avgFlow( k-1,: );
        avgFlow = [avgFlow; tmp];
    end 
end

function importantPts = ThresholdSaliency(saliencyMap)
    
    nFrames = size(saliencyMap,3);
    importantPts = [];
    for k = 1 : nFrames
        frame_map = saliencyMap(:,:,k);
        frame_map = frame_map / max(max(frame_map));
        [r c] = find(frame_map>0.75);
        sizer = size(r);
        importantPts = [importantPts ; [repmat(k,sizer,1) r c]];
    end

end

% function CropWithHomography(FolderName, FileName , CROP , RATIO , mode) %mode = 'optical', 'saliency' or 'homography'
% 
%     Reader = ReadFunctions;
%     [video,nFrames,vidHeight,vidWidth,vidFPS] = Reader.ReadData(FolderName,FileName);
%     fprintf('Number of frames = %d\n',nFrames);
% 
%     % Create movie structure & find optical flow per shot
%     mov = Reader.NewMovie(nFrames , vidHeight   ,vidWidth);
%     mov = Reader.ReadMovie(mov , video );
% 
%     Video = VideoFunctions;
%     saliencyPoints = Video.CalculateStaticSaliency(mov);
%     avgSaliency = Video.CalculateMeanSaliency( nFrames , saliencyPoints );
%     shotBoundaries = Video.DetectShotBoundaries( mov );
%     [avgFlowOptical_h] = Video.CreateFlow( mode , shotBoundaries , avgSaliency , mov , 4 ); % 'saliency' or 'optical' or 'homography'
% 
%     Cropper = CropFunctions;
%     cropArray = Cropper.EstimateWindowSize(avgFlowOptical_h, saliencyPoints , shotBoundaries , CROP , RATIO , [vidWidth vidHeight]);
%     [cropX, cropY] = Cropper.CreateWindow( cropArray .* RATIO(1,1) , cropArray .* RATIO(1,2) , avgFlowOptical_h(1:nFrames,:) , vidWidth , vidHeight );
%     
%     % Create cropped movie with cropX cropY values
%     movCropped = Reader.NewMovie(nFrames , vidHeight   ,vidWidth);
%     for k = 1 : nFrames
%         crop = mov(k).cdata( cropY(k,1):cropY(k,2) , cropX(k,1):cropX(k,2) , : );
%         movCropped(k).cdata = imresize(crop , [CROP.*RATIO(1,2)*2 ,CROP.*RATIO(1,1)*2]);
%     end
% 
%     % Save the output
%     movie2avi(movCropped, strcat(FolderName , FileName,'_cropped_optical.avi') , 'fps',  vidFPS);
% 
% end
