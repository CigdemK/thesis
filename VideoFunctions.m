function F = VideoFunctions
    F.CalculateMeanSaliency = @CalculateMeanSaliency;
    F.DetectShotBoundaries = @DetectShotBoundaries;
    F.CreateFlow = @CreateFlow;
    F.CalculateStaticSaliency = @CalculateStaticSaliency;
    F.CalculateVideoSaliency = @CalculateVideoSaliency;
    F.OpticalFlowMap = @OpticalFlowMap;
    F.CalculateHomography = @CalculateHomography;
end


%% PUBLIC FUNCTIONS

function [avg] = CalculateMeanSaliency( maxFrame  , eyes )
    avg = [];
    for k = 1 : maxFrame
        indices = find( eyes( : , 1 ) == k );
        avg = [ avg ; mean( eyes ( indices(:) , 2:3 ) ) ];
    end   
end

function [shotBoundaries] = DetectShotBoundaries( mov ) 

    nrFrames = length(mov);
    previous = imhist(mov(1,1).cdata(:));
    
    for k = 2 : nrFrames
        current = imhist(mov(1,k).cdata(:));
        diff = abs(current-previous);
        absDiff(k) = sum(diff(1));
        previous = current;
    end
    
    Util = UtilFunctions;
    shotBoundaries = Util.peakfinder(absDiff , 6*mean(absDiff));
    
    if isempty(shotBoundaries)
        shotBoundaries = [1;nrFrames];
    elseif shotBoundaries(size(shotBoundaries)) ~= nrFrames
        shotBoundaries = [1;shotBoundaries';nrFrames];
    else
        shotBoundaries = [1;shotBoundaries'];
    end
end

function [avgFlow] = CreateFlow(mode, shotBoundaries, avgSaliency , mov , WINDOW_SIZE , opticalFlowMap )

    if nargin < 6
        [opticalFlowX,opticalFlowY] = OpticalFlowMap(mov);
        opticalFlowMap(:,:,:,1) = opticalFlowX;
        opticalFlowMap(:,:,:,2) = opticalFlowY;
    end
    
    avgFlow = [];
    nrShots = length(shotBoundaries)-1;
    
    for k = 1:nrShots
        shotStart = shotBoundaries(k);
        shotEnd = shotBoundaries(k+1)-1;

        shotSaliency = avgSaliency( shotStart:shotEnd , :);
        
        if strcmp(mode,'optical')
            [Vx,Vy] = CameraSpeed(mov,opticalFlowMap);
            avgFlowTmp = CalculateFlowMean( shotSaliency(1,:) , Vx ,Vy , 'OpticalFlow' );   
        elseif strcmp(mode,'saliency')  
            [Vx,Vy] = SaliencySpeed(shotSaliency);
            avgFlowTmp = CalculateFlowMean( shotSaliency(1,:) , Vx ,Vy , 'OpticalFlow' );   
        elseif strcmp(mode,'homography')
            classes = CameraMotion(mov(shotStart:shotEnd) , WINDOW_SIZE );
            avgFlowTmp = CalculateFlowMean(shotSaliency(1,:) , 0 , 0 , 'CameraMotion' , classes );
        end 

        avgFlow = [avgFlow;avgFlowTmp];
    end

end

function [saliencyPoints , saliencyMap] = CalculateStaticSaliency(mov,algorithm)
    
    if nargin < 2
        algorithm = 'Judd';
    end
    
    saliencyPoints = [];
    nFrames = length(mov);
    
    saliencyMap = StaticSaliency(mov,algorithm);
    
    for k = 1 : nFrames
        frame_map = saliencyMap(:,:,k);
        [r c] = find(frame_map>0);
        sizer = size(r);
        saliencyPoints = [saliencyPoints ; [repmat(k,sizer,1) r c]];

    end
%     save('../CAMO_Videos/Zoom/Zoom0011.avi.mat','saliencyPoints','-append');
end

function [saliencyPoints, opticalFlowMap ] = CalculateVideoSaliency(mov)
        
    nFrames = length(mov);
    saliencyPoints = [];
    [saliencyMap , opticalFlowMap ] = VideoSaliency(mov);

    saliencyMap = mat2gray(saliencyMap);
   
    for k = 1 : nFrames
        frame_map = saliencyMap(:,:,k);
        frame_map = frame_map / max(max(frame_map));
        [r c] = find(frame_map>0.7);
        sizer = size(r);
        saliencyPoints = [saliencyPoints ; [repmat(k,sizer,1) r c]];
    end
    
end

function [opticalFlowX,opticalFlowY] = OpticalFlowMap(mov)

    disp('Calculating optical flow...');
    tic;

    nFrames = length(mov);
    [vidHeight vidWidth ~]= size(mov(1).cdata);
    opticalFlowX = zeros(vidHeight,vidWidth,nFrames);
    opticalFlowY = zeros(vidHeight,vidWidth,nFrames);
    tic;
    for k = 1:nFrames-1
        [Vx  Vy ~] = Coarse2FineTwoFrames(mov(k).cdata,mov(k+1).cdata);
        opticalFlowX(:,:,k) = Vx; 
        opticalFlowY(:,:,k) = Vy; 
    end
    
    toc;
%     save('../CAMO_Videos/Tilt/Tilt0002.avi.mat','opticalFlowX','opticalFlowY','-append');
%     load('../CAMO_Videos/Tilt/Tilt0007.avi.mat','opticalFlowX','opticalFlowY');

end

function [H] = CalculateHomography(mov)

    disp('Calculating homography...');
    tic;

    nFrames = length(mov);
    Util = UtilFunctions;

    H = [];
    for n = 1 : nFrames-1
        im1 = mov(n).cdata;
        im2 = mov(n+1).cdata;
        try
            h = Util.Homography( im1 , im2 );
        catch err
            h = zeros(3);   
        end
        H = cat(3,H,h);
        
    end
    H = cat(3,H,h);
    toc;
    
%     save('../CAMO_Videos/Dolly/Dolly0011.avi.mat','H','-append');
%     load('../CAMO_Videos/Tilt/Tilt0007.avi.mat','H');

end
%% PRIVATE FUNCTIONS

function [Vx,Vy] = CameraSpeed(mov,opticalFlowMap)
    
%     disp('Calculating camera speed...');
%     tic;
    
    shotLength = length(mov);
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

function [avgFlow] = CalculateFlowMean(startingPoint , Vx ,Vy , mode , classes)

    shotLength = length(Vx);
    Y = startingPoint(1,1);
    X = startingPoint(1,2);
    avgFlow = [X Y];
    meanFlow = [Vx Vy];
    
    dollyMean    = [ 1.0004   -0.0002    0.1011;   -0.0001    1.0006   -0.1093;  0  0  1];
    panMean      = [ 0.9979    0.0015   -0.7882;   -0.0009    0.9980    0.5330;  0  0  1];
    pedestalMean = [ 1.0034    0.0056   -0.5592;    0.0025    1.0089    1.2321;  0  0  1];
    truckingMean = [ 1.0006    0.0095   -1.4389;   -0.0002    1.0005    0.0251;  0  0  1];
    zoomMean     = [ 1.0038   -0.0001   -1.7892;    0.0001    1.0035   -0.4009;  0  0  1];
    tiltMean     = [ 1.0056    0.0064   -1.8622;    0.0004    1.0068    3.9544;  0  0  1];
    
    for k = 2 : shotLength %dismiss the first frame of the shot since it is set to startingPoint
        if  ( strcmp(mode, 'CameraMotion') )
            if classes(k) == 1
                CM = dollyMean;
            elseif classes(k) == 2
                CM = panMean;
            elseif classes(k) == 3
                CM = pedestalMean;
            elseif classes(k) == 4
                CM = truckingMean;
            elseif classes(k) == 5
                CM = zoomMean;
            elseif classes(k) == 6
                CM = tiltMean;
            end
        elseif ( strcmp(mode, 'OpticalFlow') )
            CM = eye(3);
        end
        transformed = CM*[meanFlow(k,:) 1]';
        tmp = transformed(1:2)' + avgFlow( k-1,: );
        avgFlow = [avgFlow; tmp];
    end 
end

function [Vx,Vy] = SaliencySpeed(avgSaliency)
    
    disp('Calculating saliency speed...');
    tic;
    
    lookahead = 4 ;
    lookback = 2 ;
    polinput =  [ 1:lookback+lookahead+1 ]' ;
    nrFrame = length(avgSaliency);
    Vx = zeros( nrFrame,1);
    Vy = zeros( nrFrame,1);
    for i = 1 : nrFrame;       
        if ( i-lookback > 0 && i+lookahead <= nrFrame)
              pX = polyfit(polinput, avgSaliency(i-lookback:i+lookahead,1 ),1) ; % X 
              pY = polyfit(polinput, avgSaliency(i-lookback:i+lookahead,2 ),1) ; % Y
              Vx(i) = pX(1) ;
              Vy(i) = pY(1) ;
        end      
    end 
    
    toc;
%     save('../CAMO_Videos/Zoom/Zoom0011.avi.mat','Vx','Vy','-append');

end

function [classes] = CameraMotion(mov , WINDOW_SIZE)

    disp('Calculating camera motion...');
    tic;
    
    load('data/f_w.mat');
    shotSize = length(mov);
    
    % Extract homography matrices of the shot
    H = CalculateHomography(mov);

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

    toc;
    % Add smoothing here
    
end

function [staticSaliencyMap] = StaticSaliency(mov, algorithm)

    disp('Calculating static saliency points...');
    tic;
    
    if nargin < 2
        algorithm = 'Judd';
    end

    nFrames = length(mov);
    [vidHeight vidWidth ~]= size(mov(1).cdata);
    staticSaliencyMap = zeros(vidHeight,vidWidth,nFrames);
    
    for k = 1 : nFrames
        
        if strcmp( algorithm , 'Judd' )
            frame_map = saliency(mov(k).cdata);
        elseif strcmp( algorithm , 'gbvs' )
            frame_map = gbvs(mov(k).cdata); 
            frame_map = frame_map.master_map_resized;
        elseif strcmp( algorithm , 'sr' )
%             img = imresize(mov(k).cdata, 64/size(mov(k).cdata, 2));
            img = mov(k).cdata;
            myFFT = fft2(img); 
            myLogAmplitude = log(abs(myFFT));
            myPhase = angle(myFFT);
            mySpectralResidual = myLogAmplitude - imfilter(myLogAmplitude, fspecial('average', 3), 'replicate'); 
            frame_map = abs(ifft2(exp(mySpectralResidual + 1i*myPhase))).^2;
            frame_map = mat2gray(imfilter(frame_map, fspecial('gaussian', [10, 10], 2.5)));
            frame_map = frame_map(:,:,1);
        end
            
        staticSaliencyMap (:,:,k) = frame_map;
        
    end
    toc;
%     save('../CAMO_Videos/Tilt/Tilt0002.avi.mat','staticSaliencyMap','-append');
%     load('../CAMO_Videos/Tilt/Tilt0007.avi.mat','staticSaliencyMap');

end

function [videoSaliencyMap , opticalFlowMap ] = VideoSaliency(mov)

    disp('Calculating video saliency points with Nyugen...');
    tic;
    
%     [vidHeight vidWidth ~]= size(mov(1).cdata);
%      nFrames = length(mov);
%     videoSaliencyMap = zeros(vidHeight , vidWidth , nFrames);
%     
%     H = CalculateHomography(mov); 
    [ opticalFlowMapX , opticalFlowMapY ] = OpticalFlowMap(mov);
    opticalFlowMap = cat(4,opticalFlowMapX,opticalFlowMapY);
%     staticSaliency = StaticSaliency(mov);
%     
%     load('data/f_w.mat');
% 
%     for n = 1 : nFrames-1
%         
%         dynamicSaliency = sqrt( opticalFlowMapX(:,:,n) .^ 2 + opticalFlowMapY(:,:,n) .^ 2 );
%         dynamicSaliency = dynamicSaliency/max(max(dynamicSaliency));
%         currentStaticSaliency = staticSaliency(:,:,n)/max(max(staticSaliency(:,:,n)));
% 
%         Divide into patches
%         i=1;
%         nrPatches = floor( vidHeight / 9 )*floor( vidWidth / 9 );
%         xy = zeros(nrPatches,2);
%         currentFrame=zeros(vidHeight,vidWidth);
%         
%         for k = 1 : floor( vidHeight / 9 )
%             for t = 1 : floor( vidWidth / 9 )
%                 xy(i,1) = k;
%                 xy(i,2) = t;
% 
%                 one   = H(1,1,n);
%                 two   = H(1,2,n);
%                 three = H(1,3,n);
%                 four  = H(2,1,n);
%                 five  = H(2,2,n);
%                 six   = H(2,3,n);
%                 seven = H(3,1,n);
%                 eight = H(3,2,n);
%                 nine  = H(3,3,n);
%                 inputs = [ one(:)   , two(:)   , three(:) , ...
%                            four(:)  , five(:)  , six(:)   , ...
%                            seven(:) , eight(:) , nine(:)  , ...
%                            xy(i,1)  , xy(i,2)]'; % input vector (6-dimensional pattern)
%                        
%                 wi = fi_(inputs);
%                 wv = fv(inputs);
%                 wi = reshape(wi,[9 9]);
%                 wv = reshape(wv,[9 9]);
%                 
%                 wi = wi/(2*max(max(wi)));
%                 wv = wv/(2*max(max(wv)));
% 
%                 videoSaliency = wi*currentStaticSaliency((k*9-8) : (k*9) , (t*9-8) : (t*9))  + ...
%                     wv*dynamicSaliency((k*9-8) : (k*9) , (t*9-8) : (t*9) );
%                 currentFrame(k*9-8:k*9,t*9-8:t*9) = videoSaliency(:,:);
% 
%                 i = i + 1;
% 
%             end
%         end
%         
%         videoSaliencyMap(:,:,n) = currentFrame;
%     end
    load('saliencyMap.mat'); videoSaliencyMap = saliencyMap;
    toc;
%     save('../CAMO_Videos/Tilt/Tilt0002.avi.mat','videoSaliencyMap','-append');
%     load('../CAMO_Videos/Tilt/Tilt0007.avi.mat','videoSaliencyMap','opticalFlowX','opticalFlowY');

end

