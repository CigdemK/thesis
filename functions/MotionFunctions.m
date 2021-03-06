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
        H(:,:,n) = HomographyHelper( frames(:,:,:,n), frames(:,:,:,n+1));
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
        
        VxTemp = mean(mean([(VxTemp(:)<meanX*0.7) (VxTemp(:)>meanX*(-0.7))]));
        VyTemp = mean(mean([(VyTemp(:)<meanY*0.7) (VyTemp(:)>meanY*(-0.7))]));
        Vx(t) = VxTemp;    
        Vy(t) = VyTemp;   

    end
    Vx(t+1) = Vx(t);
    Vy(t+1) = Vy(t);
       
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
% Taken from VLFeat website


    % Apply SIFT and find matching points
    [f1,d1] = vl_sift(single(rgb2gray(img1))) ;
    [f2,d2] = vl_sift(single(rgb2gray(img2))) ;

    [matches, ~] = vl_ubcmatch(d1,d2) ;

    numMatches = size(matches,2) ;

    X1 = f1(1:2,matches(1,:)) ; X1(3,:) = 1 ;
    X2 = f2(1:2,matches(2,:)) ; X2(3,:) = 1 ;

    % Find transformation matrix and apply transformation
    clear H score ok ;
    for t = 1:100
        % estimate homograpyh
        subset = vl_colsubset(1:numMatches, 4) ;
        A = [] ;
        for i = subset
            A = cat(1, A, kron(X1(:,i)', vl_hat(X2(:,i)))) ;
        end
        [~,~,V] = svd(A) ;
        H{t} = reshape(V(:,9),3,3) ;

        % score homography
        X2_ = H{t} * X1 ;
        du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:) ;
        dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:) ;
        ok{t} = (du.*du + dv.*dv) < 6*6 ;
        score(t) = sum(ok{t}) ;
    end

    [~, best] = max(score) ;
    H = H{best} ;

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
