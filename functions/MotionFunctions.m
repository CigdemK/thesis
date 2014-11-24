function F = MotionFunctions
    F.MyHomography = @MyHomography;
    F.MyOpticalFlow = @MyOpticalFlow;
    F.SaveHomographyData = @SaveHomographyData;
    F.SaveOpticalFlowData = @SaveOpticalFlowData;   
    F.CameraSpeed = @CameraSpeed;
    F.CameraMotionClass = @CameraMotionClass;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Public Functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H = MyHomography(frames)

    nFrames = size(frames,4);
    H = zeros(3,3,nFrames);
    
    for n = 1 : nFrames-1
        try
            H(:,:,n) = HomographyHelper( frames(:,:,:,n), frames(:,:,:,n+1));
        catch err
            H(:,:,n) = zeros(3);   
        end
    end
    H(:,:,end) = H(:,:,end-1);
    
end

function [opticalFlowX,opticalFlowY] = MyOpticalFlow(frames)

    [vidHeight, vidWidth, ~, nFrames]= size(frames);
    
    opticalFlowX = zeros(vidHeight,vidWidth,nFrames);
    opticalFlowY = zeros(vidHeight,vidWidth,nFrames);

    para = OpticalParameterPrep();
    for k = 1:nFrames-1
        [Vx, Vy, ~] = Coarse2FineTwoFrames(frames(:,:,:,k),frames(:,:,:,k+1),para);
        opticalFlowX(:,:,k) = Vx; 
        opticalFlowY(:,:,k) = Vy; 
    end
  
    opticalFlowX(:,:,end) = opticalFlowX(:,:,end-1);
    opticalFlowY(:,:,end) = opticalFlowY(:,:,end-1);
   
end

function SaveOpticalFlowData(folderName, movieNames)

    for i = 1:size(movieNames,2)
        
        fileName = movieNames{i};

        video = VideoReader([folderName fileName]);
        frames = read(video);
        
        [flowX,flowY] = MyOpticalFlow( frames ); 
        opt = cat(4,flowX,flowY);
        save([folderName fileName '_OpticalFlow.mat'],'opt');
        
    end
    
end

function SaveHomographyData(folderName, movieNames)
       
    for i = 1:size(movieNames,2)
        
        fileName = movieNames{i};

        video = VideoReader([folderName fileName]);
        frames = read(video);
        
        h = Video.CalculateHomography( frames ); 
        save([folderName fileName '_HomographyMatrix.mat'],'h');
        
    end
    
end

function [Vx,Vy] = CameraSpeed(opticalFlowMap)
    
%     disp('Calculating camera speed...');
%     tic;
    
    shotLength = size(opticalFlowMap,3);
    Vx = zeros(shotLength,1);
    Vy = zeros(shotLength,1);

    opticalFlowX = opticalFlowMap(:,:,:,1);
    opticalFlowY = opticalFlowMap(:,:,:,2);

    for t = 1 : shotLength-1
        
        VxTemp = opticalFlowX(:,:,t);
        VyTemp = opticalFlowY(:,:,t);
        
        meanX = mean(VxTemp(:));
        meanY = mean(VyTemp(:));
        
        VxTemp = mean(mean([(VxTemp(:)<meanX*0.9) (VxTemp(:)>meanX*(-0.9))]));
        VyTemp = mean(mean([(VyTemp(:)<meanY*0.9) (VyTemp(:)>meanY*(-0.9))]));
        Vx(t) = VxTemp;    
        Vy(t) = VyTemp;   

    end
    
    lastX = Vx( length( Vx ) );
    lastY = Vy( length( Vy ) );
    Vx(t+1) = lastX;
    Vy(t+1) = lastY;
   
%     toc;
       
end

function [classes] = CameraMotionClass(frames , WINDOW_SIZE)

%     disp('Calculating camera motion...');
%     tic;
    
    load('data/f_w.mat');
    shotSize = size(frames,4);
    
    % Extract homography matrices of the shot
    H = MyHomography(frames);

    classes = zeros(shotSize-WINDOW_SIZE,1);
    for n = 1+floor(WINDOW_SIZE/2) : shotSize-floor(WINDOW_SIZE/2)
        nnInput = [];
        for j = -WINDOW_SIZE/2:WINDOW_SIZE/2-1
            nnInput = [nnInput ; H(1,1,n+j) ; H(1,2,n+j) ; ...
                            H(1,3,n+j); H(2,1,n+j) ; H(2,2,n+j) ; ...
                            H(2,3,n+j); H(3,1,n+j) ; H(3,2,n+j)];
        end 
        
        % Check if the current camera motion is tilt
        if svmclassify(svmTilt,nnInput') == 2 
            classes(n-WINDOW_SIZE/2) = 6; 
        else
            [nnInput , ~ ] = mapminmax(nnInput,-1,+1);
            nnResult =  net(nnInput);
            classes(n-WINDOW_SIZE/2) = find(nnResult==max(nnResult));
        end
                
    end
    
    for n = 1:floor(WINDOW_SIZE/2)
        classes = [ classes( 1 ) ; classes ;classes(shotSize-WINDOW_SIZE) ];
    end

%     toc;
    % Add smoothing here
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Private Functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ H ] = HomographyHelper(img1,img2)
% Maps the points in the second image to the first image.
% Returns the homography matrix of img1-->img2

    I = single(rgb2gray(img1));
    J = single(rgb2gray(img2));

    % Apply SIFT and find matching points
    [features1,distances1] = vl_sift(I);
    [features2,distances2] = vl_sift(J);
    [matches, ~] = vl_ubcmatch(distances1,distances2,5);

    % Find transformation matrix and apply transformation
    H = findHomography( ...
         [features1( 1 , matches( 1 , : ) ) ; features1( 2 , matches( 1 , : ) ) ],...
         [features2(1 , matches( 2 , : ));features2( 2 , matches( 2 , : ) ) ]);

end

function [ params ] = OpticalParameterPrep
    alpha = 0.012;
    ratio = 0.75;
    minWidth = 20;
    nOuterFPIterations = 7;
    nInnerFPIterations = 1;
    nSORIterations = 30;
    params = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];
end
