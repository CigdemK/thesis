function F = ImageSaliency
    F.MyGbvs = @MyGbvs;
    F.MyJudd = @MyJudd;
    F.MyItti = @MyItti;
    F.MySR = @MySR;
    F.MyGaussByOpticalFlow = @MyGaussByOpticalFlow;
    F.ImpTrj = @ImpTrj;
    F.ActionsInTheEyeFixation = @ActionsInTheEyeFixation;
end

function saliencyMap = MyItti(moviePath,frames)
    if nargin < 2
        video = VideoReader( moviePath );
        frames = read(video);
    end
    
    [vidHeight, vidWidth, ~, nFrames] = size(frames);
    saliencyMap = zeros(vidHeight,vidWidth,nFrames);
    
    for i = 1:nFrames
        saliencyMap(:,:,i) = imresize(simpsal(frames(:,:,:,i)),[vidHeight vidWidth]);
    end
end

function saliencyMap = MyJudd(moviePath,frames)
       
    if nargin < 2
        video = VideoReader( moviePath );
        frames = read(video);
    end
    
    [vidHeight, vidWidth, ~, nFrames] = size(frames);
    saliencyMap = zeros(vidHeight,vidWidth,nFrames);
    
    for i = 1:nFrames
        saliencyMap(:,:,i) = saliency(frames(:,:,:,i));
    end
end

function saliencyMap = MyGbvs(moviePath,frames)
       
    if nargin < 2
        video = VideoReader( moviePath );
        frames = read(video);
    end
    
    [vidHeight, vidWidth, ~, nFrames] = size(frames);   
    saliencyMap = zeros(vidHeight,vidWidth,nFrames);
    
    for i = 1:video.NumberOfFrames
        frame_map = gbvs(frames(:,:,:,i));
        saliencyMap(:,:,i) = frame_map.master_map_resized;
    end
end

function saliencyMap = MyGaussByOpticalFlow(moviePath,frames)  

    if nargin < 2
        video = VideoReader( moviePath );
        frames = read(video);
    end
    
    [vidHeight, vidWidth, ~, nFrames] = size(frames);   
    saliencyMap = zeros(vidHeight,vidWidth,nFrames);
    opt = zeros(vidHeight,vidWidth,nFrames,2);
    
%     load([moviePath '_OpticalFlow.mat']);
    motion = MotionFunctions;
    [opt(:,:,:,1),opt(:,:,:,2)] = motion.MyOpticalFlow(frames);
    
    for k = 1:video.NumberOfFrames-1
        
        [gaussCenterX,gaussCenterY] = GetGaussCenter(opt(:,:,k,1),opt(:,:,k,2));
        gaussCenterX = gaussCenterX * video.Width;
        gaussCenterY = gaussCenterY * video.Height;

        offsetX = gaussCenterX - video.Width/2;
        offsetY = gaussCenterY - video.Height/2;

        currentGauss = MyGauss(video.Height,video.Width,video.Height/5,offsetX,offsetY);
        saliencyMap(:,:,k) = mat2gray(currentGauss);
    end
end

function saliencyMap = ImpTrj(moviePath,frames,impTrajectories)

    ImprovedT = ImprovedTrajectories;
    Util = UtilFunctions;
    Ret = RetargetMethods;
    
    if nargin < 3
        if exist([moviePath '_impTrajectories.mat'], 'file') == 2
            load([moviePath '_impTrajectories.mat']); 
        else
            [status,exeOutput] = ImprovedT.RunImprovedTrajectories(moviePath);
            if status; return; end
            impTrajectories = ImprovedT.AnalyzeOutput(exeOutput);
            save([moviePath '_impTrajectories.mat'], 'impTrajectories');
        end
    end
    
    if exist([moviePath '_staticSaliencyMap.mat'], 'file') == 2
        load([moviePath '_staticSaliencyMap.mat']); 
    else
        staticSaliency = MyJudd(moviePath,frames);
        save([moviePath '_staticSaliencyMap.mat'], 'staticSaliency');
    end
    
    [vidHeight,vidWidth,~,nFrames] = size(frames);
        
    % Get trajectories to be used for saliency map
    shotBoundaries = Util.ReadShotBoundaries(moviePath,nFrames);
    keyFrames = Ret.GetKeyFrames(staticSaliency,shotBoundaries); 
    seeds = ImprovedT.GetSeedTrjs(impTrajectories,keyFrames,shotBoundaries,...
        staticSaliency(:,:,keyFrames));
    impPts = round(Ret.GetPointsFromSeeds(impTrajectories,seeds)); 
    
    saliencyMap = zeros(vidHeight,vidWidth,nFrames);
    gaussianFilter = fspecial('gaussian',[100 100],20);
    for k = 1:nFrames

        indices = impPts( : , 2 ) == k ;
        if sum(indices) == 0; continue; end;
        
        currentSaliencyMap = zeros(vidHeight,vidWidth);
        for t = 1:size(indices,2)
            currentSaliencyMap(impPts(indices,3),impPts(indices,4))=1;
        end

        saliencyMap(:,:,k) = mat2gray(conv2(currentSaliencyMap,gaussianFilter,'same'));
    end
end

function saliencyMap = MySR(moviePath,frames)

    if nargin < 2
        video = VideoReader( moviePath );
        frames = read(video);
    end
    
    [vidHeight, vidWidth, ~, nFrames] = size(frames);   
    saliencyMap = zeros(vidHeight,vidWidth,nFrames);
    
    for k = 1 : video.NumberOfFrames
        myFFT = fft2(frames(:,:,:,k)); 
        myLogAmplitude = log(abs(myFFT));
        myPhase = angle(myFFT);
        mySpectralResidual = myLogAmplitude - ...
            imfilter(myLogAmplitude, fspecial('average', 3), 'replicate'); 
        frame_map = abs(ifft2(exp(mySpectralResidual + 1i*myPhase))).^2;
        frame_map = mat2gray(imfilter(frame_map, fspecial('gaussian', [10, 10], 2.5)));
        saliencyMap(:,:,k) = frame_map(:,:,1);       
    end
    
end

function saliencyMap = ActionsInTheEyeFixation(moviePath,frames,duration)
    
    if nargin < 3
        video = VideoReader( moviePath );
        frames = read(video);
        duration = video.Duration;
    end
    
    [vidHeight, vidWidth, ~, nFrames] = size(frames);   
    saliencyMap = zeros(vidHeight,vidWidth,nFrames);
    
    AITE = ActionsInTheEye;
    resultMap = AITE.ReadEyeTrackingData(moviePath);
    resultMap.vidDuration = duration;
    resultMap.nFrames = nFrames;
    saliencyPoints = AITE.CalculateMapping(resultMap);
    
    for k = 1:nFrames
        
        indices = find( saliencyPoints( : , 2 ) == k );
        nrIndices = length(indices);
        
        Fx = saliencyPoints( indices , 4 );
        Fy = saliencyPoints( indices , 3 );

        for i = 1:nrIndices
            saliencyMap(Fy(i),Fx(i),k) = 1;
        end
            
        gaussianFilter = fspecial('gaussian',[100 100],20);
        saliencyMap(:,:,k) = mat2gray(conv2(saliencyMap(:,:,k),gaussianFilter,'same'));
    end

end

