function F = RetargetMethods
    F.Rubinstein2009 = @Rubinstein2009;
    F.Yan2013_Mine = @Yan2013_Mine;
    F.Yan2013 = @Yan2013;
    F.LinearScaling = @LinearScaling;
    F.Cigdem = @Cigdem;
end

function retargettedFrames = Rubinstein2009(originalFrames,newSize)

    p.piecewiseThresh = 9e9;
    p.method = 'forward';
    p.seamFunc = @seamPath_gcut;
    p.s = 1;
    p.errFunc.name = @errL1;
    p.errFunc.weightNorm = @errWeightAdd;

    [vidHeight, ~, ~, nFrames] = size(originalFrames);
    
    retargettedFrames = zeros(vidHeight,newSize(2),3,nFrames);
    cd('../retargeting algorithms/Rubinstein2009');
    for i = 1:nFrames
        [retargettedFrames(:,:,:,i),~] = imretarget(originalFrames(:,:,:,i), [vidHeight newSize(2)], [], p); % vertical seams              
        toc;
    end

%     retargettedFrames = zeros(nFrames,newSize(2),3,vidHeight);
%     XTframes = permute(originalFrames,[4 2 3 1]);
%     cd('../retargeting algorithms/Rubinstein2009');
%     for i = 100:150
%         [retargettedFrames(:,:,:,i),~] = imretarget(XTframes(:,:,:,i), [nFrames newSize(2)], [], p); % vertical seams              
%         toc;
%     end
%     retargettedFrames = permute(retargettedFrames,[4 2 3 1]);


    cd('../../bin');

end

function retargettedFrames = LinearScaling(originalFrames,newSize)
    retargettedFrames = imresize(originalFrames , [newSize(1) newSize(2)]);    
end

function retargettedFrames = Yan2013_Mine(frames,newSize)

    NKP = 9;

    Crop = CropFunctions;    
    frames = Crop.RemoveBlackBars(frames);
    [vidHeight,vidWidth,~,nFrames] = size(frames);

    % Define matrices for optimized seam carving
    seamVector = zeros(vidHeight , nFrames);
    KPE = zeros( NKP , nFrames);
    KPX = zeros( NKP , nFrames);
    KPY = zeros( NKP , nFrames);
%     RPMAP = zeros(video.Height , video.Width);

    disp('Starting to remove seams...')
    oldFrames = frames;
    tic;
    for t=1:vidWidth-newSize(2); % t = number of seams removed
        
        retargettedFrames = zeros(vidHeight,vidWidth-t,3,nFrames);
        
        % For the first frame, remove the seam without optimization
        EM = findEnergy(im2double(oldFrames(:,:,:,1)));
        seamImage = findSeamImg(EM);
        seamVector(:,1) = findSeam(seamImage);  
        retargettedFrames(:,:,:,1) = SeamCut(oldFrames(:,:,:,1),seamVector(:,1));

        for n = 1:NKP

            start = floor((vidHeight/NKP)*(n-1))+1;
            eend  = floor((vidHeight/NKP)*(n));

            y = horzcat((start:eend)',seamVector(start:eend,1));
            energies = arrayfun( @( X )( EM( y( X , 1 ) , y( X , 2) ) ) , 1:size(y,1) );

            [val, ind] = max( energies );
            KPE( n , 1 ) = val;
            KPX( n , 1 ) = y( ind , 2);
            KPY( n , 1 ) = y( ind , 1);

        end
        
        % For the remaining frames

        for k = 2 : nFrames

            EM = findEnergy(im2double(oldFrames(:,:,:,k)));
            seamImage = findSeamImg(EM);
            seamVector(:,k) = findSeam(seamImage);

            for n = 1:NKP

                start = floor((vidHeight/NKP)*(n-1))+1;
                eend  = floor((vidHeight/NKP)*(n));

                y = horzcat((start:eend)',seamVector(start:eend,1));
                energies = arrayfun( @( X )( EM( y( X , 1 ) , y( X , 2) ) ) , 1:size(y,1) );

                [val, ind] = max( energies );
                KPE( n , k ) = val;
                KPX( n , k ) = y( ind , 2);
                KPY( n , k ) = y( ind , 1);

            end

            % calculate SAD for each KP (key point)
            RPMAP = zeros(vidHeight , vidWidth-t+1);
            for n = 1:NKP

                SR = vidHeight/NKP; % Search area constant
                MR = 3;             % Map area constant

                if( KPX(n,k-1)-MR < 1 );            MR = KPX(n,k-1)-1;                 end
                if( KPX(n,k-1)+MR > vidWidth-t+1);  MR = vidWidth - KPX(n,k-1)-t+1; end
                if( KPY(n,k-1)-MR < 1);             MR = KPY(n,k-1)-1;                 end
                if( KPY(n,k-1)+MR > vidHeight);     MR = vidHeight - KPY(n,k-1);    end  

                if( KPX(n,k)-SR-MR < 1 );           SR = KPX(n,k)-MR-1;                end
                if( KPX(n,k)+SR+MR > vidWidth-t+1); SR = vidWidth - KPX(n,k)-MR-t+1;end
                if( KPY(n,k)-SR-MR < 1);            SR = KPY(n,k)-MR-1;                end
                if( KPY(n,k)+SR+MR > vidHeight);    SR = vidHeight - KPY(n,k)-MR;   end      

                if( SR < MR ); MR=SR;end
                if(SR==0 || MR == 0); continue; end
                % Define search area boundaries from previous KP's
                MA_X = KPX(n,k-1)-MR : KPX(n,k-1)+MR;
                MA_Y = KPY(n,k-1)-MR : KPY(n,k-1)+MR;

                SA_X = KPX(n,k)- SR : KPX(n,k)+SR;
                SA_Y = KPY(n,k)- SR : KPY(n,k)+SR;

                for x = 1:2*SR+1    

                    MA_SA_X = floor(SA_X(x)-MR-x+1 : SA_X(x)+MR+x);
                    MA_SA_Y = floor(SA_Y(x)-MR-x+1 : SA_Y(x)+MR+x);

                    RPMAP(MA_SA_Y( 1 ):MA_SA_Y( 2*MR+1 ),MA_SA_X( 1 ):MA_SA_X( 2*MR+1 )) = ...
                        imabsdiff( EM( MA_SA_Y( 1 ):MA_SA_Y( 2*MR+1 ),MA_SA_X( 1 ):MA_SA_X( 2*MR+1 ) ) ,...
                        EM( MA_Y( 1 ):MA_Y( 2*MR+1 ) , MA_X( 1 ):MA_X( 2*MR+1 ) ) );
                end
            end

%             bot = 0;
%             top = 0;
            top = floor(max(max(RPMAP)));
            bot = ceil(min(min(RPMAP(RPMAP>0.2))));
            if( ~(top==0 & bot==0) )
                ratio = 100/(top-bot);
                RPMAP = (RPMAP-bot) .* ratio;
            end
            RPMAP(RPMAP<0)=1;

            % Adjust EM and re-compute the seam
            EM = EM .* RPMAP;
            seamImage = findSeamImg(EM);
            seamVector(:,k) = findSeam(seamImage);
            retargettedFrames(:,:,:,k) = SeamCut(oldFrames(:,:,:,k),seamVector(:,k));   
        end
        oldFrames = retargettedFrames;
        toc;
    end
end

function retargettedFrames = Yan2013(frames,newSize)

    param.alpha = 0.5;
    param.nkp = 10;
    param.mw = 3;
    param.mth = 0.2;
    group = 5;
    param.sal_interval = 10;
    
    [vidHeight, vidWidth, ~, nFrames] = size(frames);
    width_ratio =  newSize(2) / vidWidth;
    
    framesCell = cell(nFrames,1);
    for k = 1:nFrames; framesCell{k} = frames(:,:,:,k); end
    
    cd('../retargeting algorithms/Yan2013');
    framesCell = retarget_video(framesCell, width_ratio, param, group);
    cd('../../bin');
    
    retargettedFrames = zeros(vidHeight,newSize(2),3,nFrames);
    for k = 1:nFrames; retargettedFrames(:,:,:,k) = im2double(framesCell{k}); end
    
end

function retargettedFrames = Cigdem(frames,newSize)

    shotStart = 48;
    shotEnd = 95;

    Cropper = CropFunctions; 
    VideoSal = VideoSaliency;
    Util = UtilFunctions;
    ImpTrajectories = ImprovedTrajectories;
    
    frames = Cropper.RemoveBlackBars(frames);
    [vidHeight,vidWidth,~,nFrames] = size(frames);
    shotBoundaries = Util.DetectShotBoundaries(frames);
    CROP = [vidHeight/newSize(1) vidWidth/newSize(2)];
    
%     [videoSaliencyMap , opticalFlowMap] = VideoSal.Nguyen2013('',frames);
    load('videoSaliency.mat');
    load('opticalFlowMap');
    load('F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00001\ImprovedTrajectoryOriginalTJS.mat');
%     load('groups.mat');
    
%     importantPts = ThresholdSaliency(videoSaliencyMap);
    videoSaliency(videoSaliency<0.7) = 0;
    seedSaliency = videoSaliency(:,:,1) + videoSaliency(:,:,2) + videoSaliency(:,:,3);
    trajectories = ImpTrajectories.GroupTrajectories(trajectories,seedSaliency,[shotStart shotEnd]);
    avgSaliency = CalculateMeanSaliency(nFrames , importantPts );   
    [avgOpticalFlow] = CreateFlow(shotBoundaries, avgSaliency, frames, opticalFlowMap);

    cropArray = Cropper.EstimateCropWindowSize(avgOpticalFlow, importantPts, shotBoundaries, CROP);
    [cropX, cropY] = Cropper.CreateCropWindow( cropArray .* CROP(1,1) , cropArray .* CROP(1,2) , avgOpticalFlow(1:nFrames,:) , vidWidth , vidHeight );

    % Create cropped movie with cropX cropY values
    for k = 1 : nFrames
        crop = frames(cropY(k,1):cropY(k,2), cropX(k,1):cropX(k,2), :, k);
        frames(:,:,:,k) = imresize(crop , [CROP ,CROP]);
    end

%     [vidHeight, vidWidth, ~, nFrames] = size(frames);
%     load('videoSaliency.mat');
%     load('F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00001\ImprovedTrajectoryOriginalTJS.mat');
%     
%     startSaliencyMap = mean(videoSaliency(:,:,1:3),3);
%     startSaliencyMap(startSaliencyMap<0.6) = 0;
%     
%     trajectoryIndices = 1;
%     currentFrame = 1;
%     trajectoriesToBeTracked = [];
%     while 1
%         
%         currentTrajectory = trajectories{trajectoryIndices};
%         trajectoryLength = size(currentTrajectory.trajectory,2);
%         trajectoryEnd = currentTrajectory.frameNum;
%         trajectoryStart = trajectoryEnd - trajectoryLength + 1;
%         
%         if trajectoryStart == currentFrame || ...
%                 trajectoryStart == currentFrame+2 || ...
%                 trajectoryStart == currentFrame+3
%             trajectoriesToBeTracked = [trajectoriesToBeTracked;trajectoryIndices];
%         elseif trajectoryStart > currentFrame + 3
%             break;
%         end
%         trajectoryIndices = trajectoryIndices + 1;
%         
%     end
%     
%     meanX = ;
%     meanY = ; 
%     
%     
%     [vidHeight, vidWidth, ~, nFrames] = size(frames);
%     saliencyMaps = zeros(vidHeight,vidWidth,nFrames);
%     for i = 1:nFrames
%         saliencyMaps(:,:,i) = saliency(frames(:,:,:,i));
%     end
%     save('saliencyMaps.mat','saliencyMaps');
%     load('saliencyMaps.mat','saliencyMaps');
% 
%     vid = VideoSaliency;
%     videoSaliency = vid.Nguyen2013('',frames);
%     save('videoSaliency.mat','videoSaliency');
% 
%     load('videoSaliency.mat');
%     load('F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00001\ImprovedTrajectoryOriginalTJS.mat');
%     tmpsalmap = videoSaliency;
%     videoSaliency(videoSaliency<0.7) = 0;
%     imgSaliency = ImageSaliency;
%     saliencyMapFixations = imgSaliency.ActionsInTheEyeFixation('F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00001.avi');
%     
%     for i=1:10:nFrames; 
%         figure;
%         subplot(3,1,1); imshow(videoSaliency(:,:,i));
%         subplot(3,1,2); imshow(tmpsalmap(:,:,i));
%         subplot(3,1,3); imshow(saliencyMapFixations(:,:,i));
%     end
    
    
end






