function F = ImageSaliency
    F.MyGbvs = @MyGbvs;
    F.MyJudd = @MyJudd;
    F.MySR = @MySR;
    F.MyGaussByOpticalFlow = @MyGaussByOpticalFlow;
    F.Cigdem = @Cigdem;
    F.ActionsInTheEyeFixation = @ActionsInTheEyeFixation;
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

function saliencyMap = Cigdem(moviePath,trajectoriesByFrame)

    if nargin < 2
        tic;
        ImprovedT = ImprovedTrajectories;
        [status,exeOutput] = ImprovedT.RunImprovedTrajectories(moviePath);toc;
save('temp.mat');
% load('temp.mat')
        if status; return; end

        trajectories = ImprovedT.AnalyzeOutput(exeOutput);toc;
        trajectoriesByFrame = ImprovedT.GetTrajectoriesByFrame(trajectories);toc;

        regexres = regexp(moviePath,'.avi','split');
        mkdir(regexres{1});
        save([regexres{1} '\ImprovedTrajectoryOriginalTJS.mat'],'trajectories','trajectoriesByFrame');
    end
    
    video = VideoReader( moviePath );
    saliencyMap = zeros(video.Height,video.Width,video.NumberOfFrames);
    for k = 2:video.NumberOfFrames-1

        currentSaliencyMap = zeros(video.Height,video.Width);

        xIndices = round(trajectoriesByFrame{k}(1,:));
        yIndices = round(trajectoriesByFrame{k}(2,:));

        xIndices(xIndices<1) = 1;
        yIndices(yIndices<1) = 1;
        xIndices(xIndices>video.Width) = video.Width;
        yIndices(yIndices>video.Height) = video.Height;

        for t = 1:size(xIndices,2)
            currentSaliencyMap(yIndices(t),xIndices(t))=1;
        end

        gaussianFilter = fspecial('gaussian',[100 100],20);
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

function saliencyMap = ActionsInTheEyeFixation(moviePath)
    
    video = VideoReader( moviePath );
    frames = read(video);
    
    [vidHeight, vidWidth, ~, nFrames] = size(frames);   
    saliencyMap = zeros(vidHeight,vidWidth,nFrames);
    
    AITE = ActionsInTheEye;
    resultMap = AITE.ReadEyeTrackingData(moviePath);
    resultMap.vidDuration = video.Duration;
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

