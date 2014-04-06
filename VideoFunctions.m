function F = VideoFunctions
    F.CalculateMeanSaliency = @CalculateMeanSaliency;
    F.DetectShotBoundaries = @DetectShotBoundaries;
    F.CreateFlow = @CreateFlow;
    F.CalculateStaticSaliency = @CalculateStaticSaliency;
    F.CalculateVideoSaliency = @CalculateVideoSaliency;
    
    F.CalculateSaliency = @CalculateSaliency;
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
    
    shotBoundaries = peakfinder(absDiff , 6*mean(absDiff));
    
    if isempty(shotBoundaries)
        shotBoundaries = [1;nrFrames];
    elseif shotBoundaries(size(shotBoundaries)) ~= nrFrames
        shotBoundaries = [1;shotBoundaries';nrFrames];
    else
        shotBoundaries = [1;shotBoundaries'];
    end
end

function [avgFlow] = CreateFlow(mode, shotBoundaries, avgSaliency , mov , WINDOW_SIZE )

    avgFlow = zeros(length(mov),2);
    nrShots = length(shotBoundaries)-1;
    
    for k = 1:nrShots
        shotStart = shotBoundaries(k);
        shotEnd = shotBoundaries(k+1)-1;

        shotSaliency = avgSaliency( shotStart:shotEnd , :);
        
        if strcmp(mode,'optical')
            [Vx,Vy] = CameraSpeed(mov(shotStart:shotEnd));
            avgFlowTmp = CalculateFlowMean( shotSaliency(1,:) , Vx ,Vy , 'OpticalFlow' );   
        elseif strcmp(mode,'saliency')  
            [Vx,Vy] = SaliencySpeed(shotSaliency);
            avgFlowTmp = CalculateFlowMean( shotSaliency(1,:) , Vx ,Vy , 'OpticalFlow' );   
        elseif strcmp(mode,'homography')
            classes = CameraMotion(mov(shotStart:shotEnd) , WINDOW_SIZE );
%             [Vx,Vy] = CameraSpeed(mov(shotStart:shotEnd));
            avgFlowTmp = CalculateFlowMean(shotSaliency(1,:) , 0 , 0 , 'CameraMotion' , classes );
        end 
        avgFlow(shotStart:shotEnd,:) = avgFlowTmp;
    end
    avgFlowTmp = avgFlowTmp(length(avgFlowTmp),:);
%     avgFlow = [avgFlow; avgFlowTmp];
end

function [saliencyPoints] = CalculateStaticSaliency(mov)
    
    saliencyPoints = [];
    nFrames = length(mov);
    disp('Calculating saliency points...');tic;
    
    saliencyMap = StaticSaliency(mov);
    
    for k = 1 : nFrames
        frame_map = saliencyMap(:,:,k);
        [r c] = find(frame_map>0.9);
        sizer = size(r);
        saliencyPoints = [saliencyPoints ; [repmat(k,sizer,1) r c]];
        toc;
    end
%     save('../CAMO_Videos/Zoom/Zoom0011.avi.mat','saliencyPoints','-append');
end

function [saliencyPoints] = CalculateVideoSaliency(mov)

    tic;
        
    nFrames = length(mov);
    saliencyPoints = [];
    saliencyMap = VideoSaliency(mov);
    saliencyMap = mat2gray(saliencyMap);
    
    for k = 1 : nFrames
        frame_map = saliencyMap(:,:,k);
        [r c] = find(frame_map>0.9);
        sizer = size(r);
        saliencyPoints = [saliencyPoints ; [repmat(k,sizer,1) r c]];
        toc;
    end
    
end


%% PRIVATE FUNCTIONS

function [Vx,Vy] = CameraSpeed(mov)
    
        shotLength = length(mov);
        Vx = zeros(shotLength,1);
        Vy = zeros(shotLength,1);
        disp('Calculating camera speed...');tic;
            
        [opticalFlowX,opticalFlowY] = OpticalFlowMap(mov);
        
        for t = 1 : shotLength-1
            VxTemp = opticalFlowX(:,:,t);
            VyTemp = opticalFlowY(:,:,t);
            meanX = mean(VxTemp(:));
            meanY = mean(VyTemp(:));
            VxTemp = mean(mean([(VxTemp(:)<meanX*0.9) (VxTemp(:)>meanX*(-0.9))]));
            VyTemp = mean(mean([(VyTemp(:)<meanY*0.9) (VyTemp(:)>meanY*(-0.9))]));
            Vx(t) = VxTemp;    
            Vy(t) = VyTemp;   
            toc;
        end
        lastX = Vx( length( Vx ) );
        lastY = Vy( length( Vy ) );
        Vx(t+1) = lastX;
        Vy(t+1) = lastY;
       
end

function [avgFlow] = CalculateFlowMean(startingPoint , Vx ,Vy , mode , classes)

    shotLength = length(Vx);
    avgFlow = startingPoint;
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
%     save('../CAMO_Videos/Zoom/Zoom0011.avi.mat','Vx','Vy','-append');

end

function [classes] = CameraMotion(mov , WINDOW_SIZE)

    load('data/f_w.mat');
    shotSize = length(mov);
    disp('Calculating camera motion...');tic;
    
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

    % Add smoothing here
    
end

function varargout = peakfinder(x0, sel, thresh, extrema, include_endpoints)
    %PEAKFINDER Noise tolerant fast peak finding algorithm
    %   INPUTS:
    %       x0 - A real vector from the maxima will be found (required)
    %       sel - The amount above surrounding data for a peak to be
    %           identified (default = (max(x0)-min(x0))/4). Larger values mean
    %           the algorithm is more selective in finding peaks.
    %       thresh - A threshold value which peaks must be larger than to be
    %           maxima or smaller than to be minima.
    %       extrema - 1 if maxima are desired, -1 if minima are desired
    %           (default = maxima, 1)
    %       include_endpoints - If true the endpoints will be included as
    %           possible extrema otherwise they will not be included 
    %           (default = true)
    %   OUTPUTS:
    %       peakLoc - The indicies of the identified peaks in x0
    %       peakMag - The magnitude of the identified peaks
    %
    %   [peakLoc] = peakfinder(x0) returns the indicies of local maxima that
    %       are at least 1/4 the range of the data above surrounding data.
    %
    %   [peakLoc] = peakfinder(x0,sel) returns the indicies of local maxima
    %       that are at least sel above surrounding data.
    %
    %   [peakLoc] = peakfinder(x0,sel,thresh) returns the indicies of local 
    %       maxima that are at least sel above surrounding data and larger
    %       (smaller) than thresh if you are finding maxima (minima).
    %
    %   [peakLoc] = peakfinder(x0,sel,thresh,extrema) returns the maxima of the
    %       data if extrema > 0 and the minima of the data if extrema < 0
    %
    %   [peakLoc, peakMag] = peakfinder(x0,...) returns the indicies of the
    %       local maxima as well as the magnitudes of those maxima
    %
    %   If called with no output the identified maxima will be plotted along
    %       with the input data.
    %
    %   Note: If repeated values are found the first is identified as the peak
    %
    % Ex:
    % t = 0:.0001:10;
    % x = 12*sin(10*2*pi*t)-3*sin(.1*2*pi*t)+randn(1,numel(t));
    % x(1250:1255) = max(x);
    % peakfinder(x)
    %
    % Copyright Nathanael C. Yoder 2011 (nyoder@gmail.com)

    % Perform error checking and set defaults if not passed in
    error(nargchk(1,5,nargin,'struct'));
    error(nargoutchk(0,2,nargout,'struct'));

    s = size(x0);
    flipData =  s(1) < s(2);
    len0 = numel(x0);
    if len0 ~= s(1) && len0 ~= s(2)
        error('PEAKFINDER:Input','The input data must be a vector')
    elseif isempty(x0)
        varargout = {[],[]};
        return;
    end
    if ~isreal(x0)
        warning('PEAKFINDER:NotReal','Absolute value of data will be used')
        x0 = abs(x0);
    end

    if nargin < 2 || isempty(sel)
        sel = (max(x0)-min(x0))/4;
    elseif ~isnumeric(sel) || ~isreal(sel)
        sel = (max(x0)-min(x0))/4;
        warning('PEAKFINDER:InvalidSel',...
            'The selectivity must be a real scalar.  A selectivity of %.4g will be used',sel)
    elseif numel(sel) > 1
        warning('PEAKFINDER:InvalidSel',...
            'The selectivity must be a scalar.  The first selectivity value in the vector will be used.')
        sel = sel(1);
    end

    if nargin < 3 || isempty(thresh)
        thresh = [];
    elseif ~isnumeric(thresh) || ~isreal(thresh)
        thresh = [];
        warning('PEAKFINDER:InvalidThreshold',...
            'The threshold must be a real scalar. No threshold will be used.')
    elseif numel(thresh) > 1
        thresh = thresh(1);
        warning('PEAKFINDER:InvalidThreshold',...
            'The threshold must be a scalar.  The first threshold value in the vector will be used.')
    end

    if nargin < 4 || isempty(extrema)
        extrema = 1;
    else
        extrema = sign(extrema(1)); % Should only be 1 or -1 but make sure
        if extrema == 0
            error('PEAKFINDER:ZeroMaxima','Either 1 (for maxima) or -1 (for minima) must be input for extrema');
        end
    end

    if nargin < 5 || isempty(include_endpoints)
        include_endpoints = true;
    else
        include_endpoints = boolean(include_endpoints);
    end

    x0 = extrema*x0(:); % Make it so we are finding maxima regardless
    thresh = thresh*extrema; % Adjust threshold according to extrema.
    dx0 = diff(x0); % Find derivative
    dx0(dx0 == 0) = -eps; % This is so we find the first of repeated values
    ind = find(dx0(1:end-1).*dx0(2:end) < 0)+1; % Find where the derivative changes sign

    % Include endpoints in potential peaks and valleys is desired
    if include_endpoints
        x = [x0(1);x0(ind);x0(end)];
        ind = [1;ind;len0];
    else
        x = x0(ind);
    end

    % x only has the peaks, valleys, and possibly endpoints
    len = numel(x);
    minMag = min(x);


    if len > 2 % Function with peaks and valleys

        % Set initial parameters for loop
        tempMag = minMag;
        foundPeak = false;
        leftMin = minMag;

        if include_endpoints
            % Deal with first point a little differently since tacked it on
            % Calculate the sign of the derivative since we taked the first point
            %  on it does not neccessarily alternate like the rest.
            signDx = sign(diff(x(1:3)));
            if signDx(1) <= 0 % The first point is larger or equal to the second
                if signDx(1) == signDx(2) % Want alternating signs
                    x(2) = [];
                    ind(2) = [];
                    len = len-1;
                end
            else % First point is smaller than the second
                if signDx(1) == signDx(2) % Want alternating signs
                    x(1) = [];
                    ind(1) = [];
                    len = len-1;
                end
            end
        end

        % Skip the first point if it is smaller so we always start on a
        %   maxima
        if x(1) > x(2)
            ii = 0;
        else
            ii = 1;
        end

        % Preallocate max number of maxima
        maxPeaks = ceil(len/2);
        peakLoc = zeros(maxPeaks,1);
        peakMag = zeros(maxPeaks,1);
        cInd = 1;
        % Loop through extrema which should be peaks and then valleys
        while ii < len
            ii = ii+1; % This is a peak
            % Reset peak finding if we had a peak and the next peak is bigger
            %   than the last or the left min was small enough to reset.
            if foundPeak
                tempMag = minMag;
                foundPeak = false;
            end

            % Make sure we don't iterate past the length of our vector
            if ii == len
                break; % We assign the last point differently out of the loop
            end

            % Found new peak that was lager than temp mag and selectivity larger
            %   than the minimum to its left.
            if x(ii) > tempMag && x(ii) > leftMin + sel
                tempLoc = ii;
                tempMag = x(ii);
            end

            ii = ii+1; % Move onto the valley
            % Come down at least sel from peak
            if ~foundPeak && tempMag > sel + x(ii)
                foundPeak = true; % We have found a peak
                leftMin = x(ii);
                peakLoc(cInd) = tempLoc; % Add peak to index
                peakMag(cInd) = tempMag;
                cInd = cInd+1;
            elseif x(ii) < leftMin % New left minima
                leftMin = x(ii);
            end
        end

        % Check end point
        if x(end) > tempMag && x(end) > leftMin + sel
            peakLoc(cInd) = len;
            peakMag(cInd) = x(end);
            cInd = cInd + 1;
        elseif ~foundPeak && tempMag > minMag % Check if we still need to add the last point
            peakLoc(cInd) = tempLoc;
            peakMag(cInd) = tempMag;
            cInd = cInd + 1;
        end

        % Create output
        peakInds = ind(peakLoc(1:cInd-1));
        peakMags = peakMag(1:cInd-1);
    else % This is a monotone function where an endpoint is the only peak
        [peakMags,xInd] = max(x);
        if peakMags > minMag + sel
            peakInds = ind(xInd);
        else
            peakMags = [];
            peakInds = [];
        end
    end

    % Apply threshold value.  Since always finding maxima it will always be
    %   larger than the thresh.
    if ~isempty(thresh)
        m = peakMags>thresh;
        peakInds = peakInds(m);
        peakMags = peakMags(m);
    end

    % Rotate data if needed
    if flipData
        peakMags = peakMags.';
        peakInds = peakInds.';
    end

    % Change sign of data if was finding minima
    if extrema < 0
        peakMags = -peakMags;
        x0 = -x0;
    end

    % Plot if no output desired
    if nargout == 0
        if isempty(peakInds)
            disp('No significant peaks found')
        else
            figure;
            plot(1:len0,x0,'.-',peakInds,peakMags,'ro','linewidth',2);
        end
    else
        varargout = {peakInds,peakMags};
    end
end

function [opticalFlowX,opticalFlowY] = OpticalFlowMap(mov)
    disp('Calculating optical flow...');

%     nFrames = length(mov);
%     [vidHeight vidWidth ~]= size(mov(1).cdata);
%     opticalFlowX = zeros(vidHeight,vidWidth,nFrames);
%     opticalFlowY = zeros(vidHeight,vidWidth,nFrames);
%     
% %     set optical flow parameters (see Coarse2FineTwoFrames.m for the definition of the parameters)
%     alpha = 0.012;
%     ratio = 0.75;
%     minWidth = 20;
%     nOuterFPIterations = 7;
%     nInnerFPIterations = 1;
%     nSORIterations = 10;
%     para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];
% 
%     for k = 1:nFrames-1
%         [Vx  Vy ~] = Coarse2FineTwoFrames(mov(k).cdata,mov(k+1).cdata,para);
%         opticalFlowX(:,:,k) = Vx; 
%         opticalFlowY(:,:,k) = Vy; 
%         toc;
%     end
%     
%     save('../CAMO_Videos/Zoom/Zoom0011.avi.mat','opticalFlowX','opticalFlowY','-append');
    load('../CAMO_Videos/Tilt/Tilt0007.avi.mat','opticalFlowX','opticalFlowY');

end

function [H] = CalculateHomography(mov)

%     nFrames = length(mov);
%     
%     disp('Calculating homography...');
%     % Extract homography matrices f each frame
%     H = [];
%     for n = 1 : nFrames-1
%         im1 = mov(n).cdata;
%         im2 = mov(n+1).cdata;
%         try
%             cd('../CameraMotionCode');
%             [h ~] = homography( im1 , im2 );
%             cd('../bin');
%         catch err
%             h = zeros(3);   
%         end
%         H = cat(3,H,h);
%         toc;
%     end
%     H = cat(3,H,h);
%     
%     save('../CAMO_Videos/Tilt/Tilt0007.avi.mat','H','-append');
    load('../CAMO_Videos/Tilt/Tilt0007.avi.mat','H');
end

function [staticSaliencyMap] = StaticSaliency(mov)

%     disp('Calculating static saliency points...');
% 
%     nFrames = length(mov);
%     [vidHeight vidWidth ~]= size(mov(1).cdata);
%     staticSaliencyMap = zeros(vidHeight,vidWidth,nFrames);
%     
%     for k = 1 : nFrames
%         %     frame_map = gbvs(mov(k).cdata); 
%         %     [r c] = find(frame_map.master_map_resized>0.2);
%         %     frame_map = grayFrames(:,:,k);
%         frame_map = saliency(mov(k).cdata);
%         staticSaliencyMap (:,:,k) = frame_map;
%         toc;
%     end
%     
%     save('../CAMO_Videos/Tilt/Tilt0007.avi.mat','staticSaliencyMap','-append');
    load('../CAMO_Videos/Tilt/Tilt0007.avi.mat','staticSaliencyMap');

end

function [videoSaliencyMap] = VideoSaliency(mov)

%     [vidHeight vidWidth ~]= size(mov(1).cdata);
%      nFrames = length(mov);
%     videoSaliencyMap = zeros(vidHeight , vidWidth , nFrames);
%     
%     H = CalculateHomography(mov); 
%     [ opticalFlowMapX , opticalFlowMapY ] = OpticalFlowMap(mov);
%     staticSaliency = StaticSaliency(mov);
%     
%     load('data/f_w.mat');
%     disp('Calculating video saliency points...');
% 
%     for n = 1 : nFrames-1
%         
%         dynamicSaliency = sqrt( opticalFlowMapX(:,:,n) .^ 2 + opticalFlowMapY(:,:,n) .^ 2 );
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
%                 videoSaliency = wi*staticSaliency((k*9-8) : (k*9) , (t*9-8) : (t*9) , n)  + ...
%                     wv*dynamicSaliency((k*9-8) : (k*9) , (t*9-8) : (t*9) );
%                 currentFrame(k*9-8:k*9,t*9-8:t*9) = videoSaliency(:,:);
% 
%                 i = i + 1;
%                 toc;
%             end
%         end
%     
%         videoSaliencyMap(:,:,n) = currentFrame;
%     end
%     save('../CAMO_Videos/Zoom/Zoom0011.avi.mat','videoSaliencyMap','-append');
    load('../CAMO_Videos/Tilt/Tilt0007.avi.mat','videoSaliencyMap');

end

function CalculateSaliency() 
end
