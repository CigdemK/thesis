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
    cd('../retargeting algorithms/Rubinstein2009');tic;
    for i = 1:nFrames
        [retargettedFrames(:,:,:,i),~] = imretarget(originalFrames(:,:,:,i), [vidHeight newSize(2)], [], p); % vertical seams              
    end
    cd('../../bin');

end

function retargettedFrames = LinearScaling(originalFrames,newSize)
    retargettedFrames = imresize(originalFrames , [newSize(1) newSize(2)]);    
end

function retargettedFrames = Yan2013_Mine(frames,newSize)

    NKP = 10;

    Crop = CropFunctions;    
    frames = Crop.RemoveBlackBars(frames);
    [vidHeight,vidWidth,~,nFrames] = size(frames);

    % Define matrices for optimized seam carving
    seamVector = zeros(vidHeight , nFrames);
    KPE = zeros( NKP , nFrames);
    KPX = zeros( NKP , nFrames);
    KPY = zeros( NKP , nFrames);
%     RPMAP = zeros(video.Height , video.Width);

    saliencyMap = zeros(vidHeight,vidWidth,nFrames);
    for i = 1:nFrames; saliencyMap(:,:,i) = IttiSalMap(frames(:,:,:,i)); end;

    oldFrames = frames;
    oldSaliencyMap = saliencyMap;
    disp('Starting to remove seams...')
    for t=1:vidWidth-newSize(2); % t = number of seams removed
        
        retargettedFrames = zeros(vidHeight,vidWidth-t,3,nFrames);
        retargettedSaliency = zeros(vidHeight,vidWidth-t,nFrames);
        
        % For the first frame, remove the seam without optimization
        
        EM = FindEnergyYan(mat2gray(oldFrames(:,:,:,1)),oldSaliencyMap(:,:,i));
        seamImage = findSeamImg(EM);
        seamVector(:,1) = findSeam(seamImage);
        retargettedFrames(:,:,:,1) = SeamCut(oldFrames(:,:,:,1),seamVector(:,1));
        retargettedSaliency(:,:,1) = SeamCut(oldSaliencyMap(:,:,1),seamVector(:,1));

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

            EM = FindEnergyYan(im2double(oldFrames(:,:,:,k)),oldSaliencyMap(:,:,i));
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

                if( KPX(n,k-1)-MR < 1 );            MR = KPX(n,k-1)-1;              end
                if( KPX(n,k-1)+MR > vidWidth-t+1);  MR = vidWidth - KPX(n,k-1)-t+1; end
                if( KPY(n,k-1)-MR < 1);             MR = KPY(n,k-1)-1;              end
                if( KPY(n,k-1)+MR > vidHeight);     MR = vidHeight - KPY(n,k-1);    end  

                if( KPX(n,k)-SR-MR < 1 );           SR = KPX(n,k)-MR-1;             end
                if( KPX(n,k)+SR+MR > vidWidth-t+1); SR = vidWidth - KPX(n,k)-MR-t+1;end
                if( KPY(n,k)-SR-MR < 1);            SR = KPY(n,k)-MR-1;             end
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
            retargettedSaliency(:,:,k) = SeamCut(oldSaliencyMap(:,:,k),seamVector(:,k));
        end
        oldFrames = retargettedFrames;
        oldSaliencyMap = retargettedSaliency;
    end
    retargettedFrames = uint8(retargettedFrames);
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

function retargettedFrames = Cigdem(videoPath,frames,newSize)

    % Init Definitions
    Cropper = CropFunctions; 
    VideoSal = VideoSaliency;
    Util = UtilFunctions;
    ImpTrjs = ImprovedTrajectories;
    Ali = Ali59MFunctions;
    
    % Preprocess input video & definitions
    frames = Cropper.RemoveBlackBars(frames);
    [vidHeight,vidWidth,~,nFrames] = size(frames);
    shotBoundaries = Util.ReadShotBoundaries(videoPath,nFrames);
    cropRatio = [newSize(1)/vidHeight newSize(2)/vidWidth];
 
%     % Calculate saliency & optical flow maps & trajectories
%     [videoSaliency , opticalFlow, staticSaliency] = VideoSal.Nguyen2013('',frames);
%     [status,exeOutput] = ImpTrjs.RunImprovedTrajectories(videoPath);toc;
% 
%     if ~status
%         impTrajectories = ImpTrjs.AnalyzeOutput(exeOutput);
%     end
%         
%     save([videoPath '_videoSaliencyMap.mat'], 'videoSaliency'); 
%     save([videoPath '_opticalFlow.mat'], 'opticalFlow');
%     save([videoPath '_staticSaliencyMap.mat'], 'staticSaliency'); 
%     save([videoPath '_impTrajectories.mat'], 'impTrajectories');
    load([videoPath '_videoSaliencyMap.mat']); 
    load([videoPath '_opticalFlow.mat']);
    load([videoPath '_staticSaliencyMap.mat']); 
    load([videoPath '_impTrajectories.mat']);

    % Get important parts with key frames and trajectories
    keyFrames = GetKeyFrames(staticSaliency,shotBoundaries); 

    seeds     = ImpTrjs.GetSeedTrjs(impTrajectories,keyFrames,shotBoundaries,staticSaliency(:,:,keyFrames)); toc;
    impPts    = GetPointsFromSeeds(impTrajectories,seeds); toc;

    % Apply cropping
    avgSaliency     = Ali.CalculateMeanSaliency(impPts, shotBoundaries );  toc; 
    avgOpticalFlow  = Ali.CreateFlow(shotBoundaries, avgSaliency, frames, opticalFlow); toc;
    cropShots       = Cropper.GetWindowSize(avgOpticalFlow, impPts, shotBoundaries, newSize, [vidWidth vidHeight]);  toc;
    [cropX,cropY]   = Cropper.CreateWindow(cropShots, shotBoundaries, avgOpticalFlow, [vidWidth vidHeight]); toc;
    
    [~,maxShotSize]=(max(cropShots(:,1)));
    newSize = [cropShots(maxShotSize,2) cropShots(maxShotSize,1)];
    % Create cropped movie with cropX cropY values
    for k = 1 : nFrames
        retargettedFrames(:,:,:,k) = imresize(frames(cropY(k,1):cropY(k,2),...
            cropX(k,1):cropX(k,2), :, k),newSize);
    end

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Private Functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pts = GetPointsFromSeeds(trajectories,seeds)

    nSeeds = size(seeds,1);
    
    seedLengths = arrayfun(@(X)size(trajectories(seeds(X)).trajectory,2),1:nSeeds);
    seedEnds = [trajectories(seeds).frameNum];
    seedStarts = seedEnds - seedLengths + 1;
    
    pts = [];
    for t = 1:nSeeds
        
        currentTrajectory = trajectories(seeds(t)).trajectory';
        currentStart = seedStarts(t);
        currentLength = seedLengths(t);
        
        pts = [pts; [repmat([0 currentStart],[currentLength,1]) currentTrajectory]];
            
    end
end

function keyFrames = GetKeyFrames(saliencyMaps,shotBoundaries)
    
    nrOfShots = size(shotBoundaries,1) -1;
    keyFrames = [];

    for i = 1:nrOfShots
        shotSize = shotBoundaries(i+1)-shotBoundaries(i)+1;
        currentShotSaliency = saliencyMaps(:,:,shotBoundaries(i):shotBoundaries(i+1));
        dispersion = CalculateDispersion(currentShotSaliency);
        peakFrames = [];
        orderOfFit = 3;
        while isempty(peakFrames)
            dispersionPolyCoeff = polyfit(1:shotSize,dispersion',orderOfFit);
            dispersionPolyVal = polyval(dispersionPolyCoeff,1:shotSize);
            [~,peakFrames] = findpeaks(dispersionPolyVal);
            orderOfFit = orderOfFit + 1;
        end
        keyFrames = [keyFrames;peakFrames(:)+shotBoundaries(i)-1];
    end

end

function dispersion = CalculateDispersion(videoSaliency,threshold)

    if nargin < 2
        threshold = 0.8;
    end
    
    nFrames = size(videoSaliency,3);
    saliencyMaps = mat2gray(videoSaliency);
    dispersion = zeros(nFrames,1);
    salPts = [];
    for i = 1:nFrames 

        currentSaliency = saliencyMaps(:,:,i);
        [x,y] = find(currentSaliency > threshold);
        nrOfPoints = size(x,1);
    %     X(i-shotStart+1) = var(x);
    
        currentDispersion = 0;
        for k = 1: nrOfPoints
            for t = k+1: nrOfPoints
                currentDispersion = currentDispersion + (double((x(k)-x(t))^2 + (y(k)-y(t))^2));
            end
        end
        dispersion(i) = currentDispersion / (nrOfPoints^2);
    end
end

function sal = IttiSalMap(img)
    [imHeight, imWidth, ~] = size(img);
    sal = simpsal(img);
    sal = imresize(sal, [imHeight imWidth]);
end

function EM = FindEnergyYan(x,sal)

    [imHeight, imWidth, imDim] = size(x);

    Grdx =[ 1   2   1;
            0   0   0;
           -1  -2  -1];
    Grdy =[-1   0   1;
           -2   0   2;
           -1   0   1];

    for i=1:imDim
        Eh(:,:,i)=conv2(x(:,:,i),Grdx,'same');
        Ev(:,:,i)=conv2(x(:,:,i),Grdy,'same');
        E(:,:,i)=sqrt(Eh(:,:,i).^2+Ev(:,:,i).^2);
    end
    gr = 1 / imDim*sum(E,3);   %finds average gradient image if RGB image
    
    EM = 0.5.*sal + 0.5.*gr ;
    
end



