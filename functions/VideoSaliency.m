function F = VideoSaliency
    F.Nguyen2013 = @Nguyen2013;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PUBLIC FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [videoSaliencyMap , opticalFlow, staticSaliency] = Nguyen2013(moviePath,frames)

    if nargin < 2
        video = VideoReader( moviePath );
        frames = read(video);
    end
    
    [vidHeight, vidWidth, ~, nFrames] = size(frames);
    videoSaliencyMap = zeros(vidHeight,vidWidth,nFrames);
    
    Motion = MotionFunctions;tic;
    H = Motion.MyHomography(frames);
    [opticalFlowMapX, opticalFlowMapY] = Motion.MyOpticalFlow(frames);
    opticalFlow = cat(4,opticalFlowMapX,opticalFlowMapY);toc;
%     load([moviePath '_opticalFlow.mat']);
%     opticalFlowMapX = opticalFlow(:,:,:,1);
%     opticalFlowMapY = opticalFlow(:,:,:,2);
    StaticSal = ImageSaliency;
    staticSaliency = StaticSal.MyJudd('',frames);toc;
%     load([moviePath '_staticSaliencyMap.mat']);
    
    load('../myData/f_w.mat');
    for n = 1 : nFrames-1
        
        dynamicSaliency = sqrt( opticalFlowMapX(:,:,n) .^ 2 + opticalFlowMapY(:,:,n) .^ 2 );
        dynamicSaliency = dynamicSaliency/max(max(dynamicSaliency));
        currentStaticSaliency = staticSaliency(:,:,n)/max(max(staticSaliency(:,:,n)));

%       Divide into patches
        i=1;
        nrPatches = floor( vidHeight / 9 )*floor( vidWidth / 9 );
        xy = zeros(nrPatches,2);
        currentFrame=zeros(vidHeight,vidWidth);
        
        for k = 1 : floor( vidHeight / 9 )
            for t = 1 : floor( vidWidth / 9 )
                xy(i,1) = k;
                xy(i,2) = t;

                one   = H(1,1,n);
                two   = H(1,2,n);
                three = H(1,3,n);
                four  = H(2,1,n);
                five  = H(2,2,n);
                six   = H(2,3,n);
                seven = H(3,1,n);
                eight = H(3,2,n);
                nine  = H(3,3,n);
                inputs = [ one(:)   , two(:)   , three(:) , ...
                           four(:)  , five(:)  , six(:)   , ...
                           seven(:) , eight(:) , nine(:)  , ...
                           xy(i,1)  , xy(i,2)]'; % input vector (6-dimensional pattern)
                       
                wi = fi_(inputs);
                wv = fv(inputs);
                wi = reshape(wi,[9 9]);
                wv = reshape(wv,[9 9]);
                
                wi = wi/(2*max(max(wi)));
                wv = wv/(2*max(max(wv)));

                videoSaliency = wi*currentStaticSaliency((k*9-8) : (k*9) , (t*9-8) : (t*9))  + ...
                    wv*dynamicSaliency((k*9-8) : (k*9) , (t*9-8) : (t*9) );
                currentFrame(k*9-8:k*9,t*9-8:t*9) = videoSaliency(:,:);

                i = i + 1;

            end
        end
        videoSaliencyMap(:,:,n) = currentFrame;
    end
    videoSaliencyMap = mat2gray(videoSaliencyMap); toc;
    
end


% function [saliencyPoints, opticalFlowMap ] = Nguyen2013(moviePath)
% 
%     saliencyPoints = [];
%     [saliencyMap , opticalFlowMap ] = Nguyen2013Core(moviePath);
% 
%     saliencyMap = mat2gray(saliencyMap);
%    
%     
%     for k = 1 : nFrames
%         frame_map = saliencyMap(:,:,k);
%         frame_map = frame_map / max(max(frame_map));
%         [r, c] = find(frame_map>0.75);
%         sizer = size(r);
%         saliencyPoints = [saliencyPoints ; [repmat(k,sizer,1) r c]];
%     end
%     
% end
