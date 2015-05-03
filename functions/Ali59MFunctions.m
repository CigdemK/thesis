function F = Ali59MFunctions
    F.ShowSaliencyPoints = @ShowSaliencyPoints;
    F.CropWithVideoSaliency = @CropWithVideoSaliency;
    F.CalculateMeanSaliency = @CalculateMeanSaliency;
    F.CreateFlow = @CreateFlow;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Public Functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function frames = ShowSaliencyPoints(moviePath , mode, saliency, saliencyPoints) 
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

    Cropper = CropFunctions;
    if strcmp( mode , 'all' )
        
        for k = 1 : video.NumberOfFrames

            indices = find( saliencyPoints( : , 2 ) == k );
            nrIndices = size(indices,1);
            
            if nrIndices == 0; continue; end;

            Fx = saliencyPoints( indices , 4 );
            Fy = saliencyPoints( indices , 3 );
            redDot = cat( 3 , 254 * ones(6),zeros(6),zeros(6));
            
            for i = 1:nrIndices
                [cropX, cropY] = Cropper.CreateWindow([5 5] , [1 2], [Fx(i) Fy(i)] , [video.Width video.Height]);
                frames(cropY(1,1):cropY(1,2) ,cropX(1,1):cropX(1,2) ,:, k) = redDot ;   
            end
      
        end
        
    elseif strcmp( mode , 'mean' )
        
        if strcmp( saliency , 'given' ) 
            avgSaliency = saliencyPoints;
        else
            avgSaliency = CalculateMeanSaliency(video.NumberOfFrames ,saliencyPoints);
        end
        
        redDot = cat( 3 , 254 * ones(16),zeros(16),zeros(16));
        
        for k = 1 : video.NumberOfFrames
            
            Fx = avgSaliency(k,2);
            Fy = avgSaliency(k,1);
            [cropX, cropY] = Cropper.CreateWindow([15 15],[1 2],[Fx  Fy] , [video.Width  video.Height] );
            frames(cropY(1,1):cropY(1,2) ,cropX(1,1):cropX(1,2) ,:, k) = redDot ;
        
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

    cropArray = Cropper.GetWindowSize(avgOpticalFlow, importantPts, shotBoundaries, CROP);
    [cropX, cropY] = Cropper.CreateWindow( cropArray .* CROP(1,1) , cropArray .* CROP(1,2) , avgOpticalFlow(1:nFrames,:) , vidWidth , vidHeight );

    % Create cropped movie with cropX cropY values
    for k = 1 : nFrames
        crop = frames(cropY(k,1):cropY(k,2), cropX(k,1):cropX(k,2), :, k);
        frames(:,:,:,k) = imresize(crop , [CROP ,CROP]);
    end

end

function avgKeys = CalculateMeanSaliency(eyes, shotBoundaries)
    avgAll = [];
    nFrames = eyes(end,2);
    nShots = size(shotBoundaries,1)-1;
    
    for k = 1 : nFrames
        indices = find( eyes( : , 2 ) == k );
        
        if isempty(indices); continue; end;
        
        avgAll = [ avgAll ; k mean( eyes ( indices(:) , 3:4 ) ) ];
    end   
    
    avgKeys = [];
    for i = 1:nShots
        availableFrame = min(avgAll(avgAll(:,1) >= shotBoundaries(i)));
        if isempty(availableFrame); continue; end;
        maskVector = avgAll(:,1)<=(availableFrame+3) & ...
            avgAll(:,1)>=(availableFrame);
        avgKeys = [avgKeys;mean(avgAll(maskVector,2:3),1)];
    end
   
end

function [avgFlow] = CreateFlow(shotBoundaries, avgSaliency , frames , opticalFlowMap )

    if nargin < 4
        Motion = MotionFunctions;
        [opticalFlowX,opticalFlowY] = Motion.MyOpticalFlow(frames);
        opticalFlowMap(:,:,:,1) = opticalFlowX;
        opticalFlowMap(:,:,:,2) = opticalFlowY;
    end
    
    avgFlow = [];
    nrShots = size(shotBoundaries,1)-1;
    [vidHeight,vidWidth,~,~] = size(frames);
    
    Motion = MotionFunctions;
    for k = 1:nrShots

        [Vx,Vy] = Motion.CameraSpeed(opticalFlowMap(:,:,shotBoundaries(k):shotBoundaries(k+1)-1,:));
        avgFlowTmp = CalculateFlowMean( avgSaliency(k,:) , Vx ,Vy);   
        avgFlowTmp(avgFlowTmp(:,1)>vidWidth,1) = vidWidth;
        avgFlowTmp(avgFlowTmp(:,2)>vidHeight,2) = vidHeight;
        avgFlowTmp(avgFlowTmp<1) = 1;
        avgFlow = [avgFlow;avgFlowTmp];
    end
    
    avgFlow = [avgFlow;avgFlow(end,:)];
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Private Functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [avgFlow] = CalculateFlowMean(startingPoint , Vx ,Vy)

    shotLength = length(Vx);
    Y = startingPoint(1,1);
    X = startingPoint(1,2);
    avgFlow = [X Y];
    meanFlow = [Vx Vy];
    
    for k = 2 : shotLength %dismiss the first frame of the shot since it is set to startingPoint
        tmp = meanFlow(k,:) + avgFlow( k-1,: );
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


