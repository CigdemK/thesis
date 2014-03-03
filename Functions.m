function F = Functions
    F.ReadData = @ReadData;
    F.ReadEyeTrackingData = @ReadEyeTrackingData;
    F.ReadMovie = @ReadMovie;
    F.NewMovie = @NewMovie;
    F.CalculateMapping = @CalculateMapping;
    F.CalculateMeanSaliency = @CalculateMeanSaliency;
    F.CreateWindow = @CreateWindow;
    F.PlotTrajectory2D = @PlotTrajectory2D;
    F.PlotTrajectory3D = @PlotTrajectory3D;
    F.DetectShotBoundaries = @DetectShotBoundaries;
    F.OpticalFlow = @OpticalFlow;
    F.CameraSpeed = @CameraSpeed;
    F.CalculateFlowMean = @CalculateFlowMean;
    F.CreateFlow = @CreateFlow;
    F.SaliencySpeed = @SaliencySpeed;
    F.FlowSpeed = @FlowSpeed;
    F.EstimateWindowSize = @EstimateWindowSize;
    F.CheckIsInside = @CheckIsInside;
    F.ChangeCropSize = @ChangeCropSize;
end

function ReadData( foldername,filename )

    % Read video
    video       = VideoReader( strcat(foldername, filename ) );
    nFrames     = video.NumberOfFrames;
    vidHeight   = video.Height;
    vidWidth    = video.Width;
    vidFPS      = video.FrameRate;
    vidDuration = video.Duration;

    save('Data.mat');   
end

function ReadEyeTrackingData(foldername,filename)

    load('Data.mat');
    % Read eye tracking data
    fid = fopen( strcat(foldername,'006_' , filename , '.txt') );
    data = textscan( fid , '%d %f %f %f %f %s' );
    fclose( fid );
    timeStamp   = cell2mat( data( 1 ) );
    xScreen     = cell2mat( data( 4 ) );
    yScreen     = cell2mat( data( 5 ) );
    
    % Read calibration of the experiment setup
    fid = fopen( '../holywood2/geometry.txt' );
    data = textscan( fid , '%f %f %f %d %d' );
    screenResX   = cell2mat( data( 4 ) );
    screenResY   = cell2mat( data( 5 ) );
    fclose( fid );
    
    % Read screen resolution data
    fid = fopen( '../holywood2/resolution.txt' );
    data = textscan( fid , '%s %d %d %f' );
    fclose( fid );
    vidNames =  [ data{ 1 } ] ;
    index = (strcmp( filename , vidNames ) );
    videoResX   = cell2mat( data( 2 ) );
    videoResX   = videoResX( find( index == 1 ) );
    videoResY   = cell2mat( data( 3 ) );
    videoResY   = videoResY( find( index == 1 ) );
   
    save('Data.mat');  
end

function [mov] = ReadMovie( mov , video , nFrames )
    for k = 1 : nFrames
        mov(k).cdata = read(video, k);    
    end
end

function [mov] = NewMovie(nFrames , vidHeight   ,vidWidth )
    mov( 1 : floor( nFrames ))= struct('cdata',zeros(vidHeight   ,vidWidth   , 3,'uint8'),'colormap',[]);
end

function CalculateMapping()
    load('Data.mat')
    a           = double( videoResX ) / double( screenResX ) ; 
    b           = ( double( screenResY ) - double( videoResY ) / a ) / 2;
    xHeight     = a * xScreen;
    yHeight     = a * ( yScreen - b );
    frames      = round( nFrames * timeStamp / ( vidDuration * 1000000 ) ) +1;
    eyes        = [ frames xHeight yHeight ];
    eyes      = eyes( eyes(:,1) <= nFrames,:  );
    save('Data.mat');
end

function [avg,distances] = CalculateMeanSaliency( maxFrame  , eyes )
    avg = [];
    distances = [];
    for k = 1 : maxFrame+1
        indices = find( eyes( : , 1 ) == k );
        avg = [ avg ; mean( eyes ( indices(:) , 2:3 ) ) ];
        if k > 1
            distances(k-1) = dist(avg(k),avg(k-1)');
        end
    end   
end

function [croppedX,croppedY] = CreateWindow( limX , limY , avg , vidWidth , vidHeight)
    
    croppedX = [];
    croppedY = [];
    for k = 1 : size( avg )
        minX = floor( avg( k , 1)) - limX(k) ;
        maxX = floor( avg( k , 1)) + limX(k) ;
        minY = floor( avg( k , 2)) - limY(k) ;
        maxY = floor( avg( k , 2)) + limY(k) ;
        if minX <= 0 
            minX = 1;
            maxX = 2 * limX(k) + 1 ;
        end
        if maxX > vidWidth 
            maxX = vidWidth ;
            minX = vidWidth - ( 2 * limX(k));
        end
        if minY <= 0 
            minY = 1;
            maxY = 2 * limY(k) + 1;
        end
        if maxY > vidHeight 
            maxY = vidHeight ;
            minY = vidHeight - ( 2 * limY(k));
        end
        croppedX = [croppedX ; [minX maxX] ];
        croppedY = [croppedY ; [minY maxY] ];
    end
    
end

function PlotTrajectory2D( data , fitPower )
    figure;
    x = double(data( : , 1 ));
    y = double(data( : , 2 ));
    p = polyfit( x , y , fitPower );
    plotx= linspace( min(x) , max(x));
    hold on;
    plot( x , y ,  'x' );    
    hold on;
    ploty = polyval( p , plotx );
    plot(plotx, ploty, '-');
    axis([0 500 0 400]);
    hold off;
end

function PlotTrajectory3D( data )
    figure;
    x = data( : , 1 );
    y = data( : , 2 );
    sizeD = length(data);
    plot3( x , y , 1:sizeD); % do a 3d trajectory mapping
    hold on;
    plot3( x , y , 1:sizeD, 'x' ); 
    hold off;
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

function [Vx,Vy] = OpticalFlow(images,alpha,iterations)
    %// Calculating optical flow of a sequence of images. 
    %// images : 3D array that contains a sequence of images. size of images is (imageHeight, imageWidth, frame number)
    %// alpha
    %// iterations. 
    [height,width,frames]=size(images);

    %//initialzation of u and v
    Vx = zeros(height,width);
    Vy = zeros(height,width);

    for k = 1:frames-1

        % //initialization of Ex Ey and Et 
        Ex = zeros(height-1,width-1,frames-1);
        Ey = zeros(height-1,width-1,frames-1);
        Et = zeros(height-1,width-1,frames-1);

        %//calculating Ex Ey and Et in frame k.
        for x = 2:width-1
            for y = 2:height-1
                Ex(y,x,k) = (images(y+1,x+1,k)-images(y+1,x,k)+images(y,x+1,k)...
                    -images(y,x,k)+images(y+1,x+1,k+1)-images(y+1,x,k+1)...
                    +images(y,x+1,k+1)-images(y,x,k+1))/4;

                Ey(y,x,k) = (images(y,x,k)-images(y+1,x,k)+images(y,x+1,k)...
                    -images(y+1,x+1,k)+images(y,x,k+1)-images(y+1,x,k+1)...
                    +images(y,x+1,k+1)-images(y+1,x+1,k+1))/4;

                Et(y,x,k) = (images(y+1,x,k+1)-images(y+1,x,k)+images(y,x,k+1)...
                    -images(y,x,k)+images(y+1,x+1,k+1)-images(y+1,x+1,k)...
                    +images(y,x+1,k+1)-images(y,x+1,k))/4;
            end
        end

        for nn = 1:iterations
            for x = 2:width-1
                for y = 2:height-1

                    Vxbar = (Vx(y-1,x)+Vx(y,x+1)+Vx(y+1,x)+Vx(y,x-1))/6+...
                         (Vx(y-1,x-1)+Vx(y-1,x+1)+Vx(y+1,x+1)+Vx(y+1,x-1))/12;

                    Vybar = (Vy(y-1,x)+Vy(y,x+1)+Vy(y+1,x)+Vy(y,x-1))/6+...
                        (Vy(y-1,x-1)+Vy(y-1,x+1)+Vy(y+1,x+1)+Vy(y+1,x-1))/12;

                    %// chapter 12 of Horn's paper
                    temp = (Ex(y,x,k)*Vxbar+Ey(y,x,k)*Vybar+Et(y,x,k))/(alpha^2 + Ex(y,x,k)^2 + Ey(y,x,k)^2);
                    %// update u and v 
                    Vx(y,x) = Vxbar-Ex(y,x,k)*temp;
                    Vy(y,x) = Vybar-Ey(y,x,k)*temp;
                end
            end
        end

    end
end

function [Vx,Vy] = CameraSpeed(mov)
    tmp = [];
    nrFrame = length(mov);
    for k = 1 : nrFrame;
        tmp2 = rgb2gray(mov(1,k).cdata);
        tmp = cat(3,tmp,double(tmp2));
    end
    [Vx , Vy] = OpticalFlow(tmp,0.001,1);
end

function [avgFlow] = CalculateFlowMean(avg , Vx ,Vy )
    time = length(Vx);
    avgFlow = avg;
    meanFlow = [Vx Vy];
    for k = 2 : time
        tmp = meanFlow(k,:) .* 1 + avgFlow( k-1,: );
        avgFlow = [avgFlow; tmp];
    end 
end

function [avgFlow] = CreateFlow(mode, shotBoundaries, avgSaliency , mov , WINDOW_SIZE )

    avgFlow = [];
    nrShots = length(shotBoundaries)-1;
    for k = 1:nrShots
        shotStart = shotBoundaries(k);
        shotEnd = shotBoundaries(k+1)-1;
        shotSaliency = avgSaliency( shotStart:shotEnd , :);
        
        if strcmp(mode,'optical')
            [Vx,Vy] = FlowSpeed(mov(shotStart:shotEnd) , WINDOW_SIZE , shotStart , shotEnd);
        elseif strcmp(mode,'saliency')  
            [Vx,Vy] = SaliencySpeed(shotSaliency);
        end
            
        avgFlowTmp = CalculateFlowMean( shotSaliency(1,:) , Vx ,Vy );        
        avgFlow = [avgFlow; avgFlowTmp];
    end
    avgFlowTmp = avgFlowTmp(length(avgFlowTmp),:);
    avgFlow = [avgFlow; avgFlowTmp];
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
end

function [Vx,Vy] = FlowSpeed(mov , WINDOW_SIZE , shotStart , shotEnd)
    
        Vx = [0];
        Vy = [0];
        shotLength = shotEnd - shotStart;
        for t = 2 : shotLength-WINDOW_SIZE 
            [VxTemp  VyTemp ] = CameraSpeed(mov(t:t+WINDOW_SIZE));
            meanX = mean(VxTemp(:));
            meanY = mean(VyTemp(:));
            VxTemp = mean(mean([(VxTemp(:)<meanX*0.9) (VxTemp(:)>meanX*(-0.9))]));
            VyTemp = mean(mean([(VyTemp(:)<meanY*0.9) (VyTemp(:)>meanY*(-0.9))]));
            Vx = [Vx; VxTemp];    
            Vy = [Vy; VyTemp];   
        end
        lastX = Vx( length( Vx ) );
        lastY = Vy( length( Vy ) );
        for t = 1:WINDOW_SIZE+1
            Vx = [Vx; lastX];
            Vy = [Vy; lastY];
        end
end

function [crop] = EstimateWindowSize(avgFlow, eyes , shotBoundaries ,minSize , ratio , maxSize)
 
    crop = [];
    nrShots = length(shotBoundaries)-1;
    for k = 1:nrShots
        shotStart = shotBoundaries(k);
        shotEnd = shotBoundaries(k+1)-1;
        shotFlow = avgFlow( shotStart:shotEnd , :);
        
        cropShot = minSize;
        for t = shotStart : shotEnd     
            
            % get all gaze pts in each frame
            indices = find( eyes( : , 1 ) == t );
            nrIndices = length(indices);
            Fx = eyes( indices( 1:nrIndices ) , 2 );
            Fy = eyes( indices( 1:nrIndices ) , 3 );
            points = [Fx Fy];
            points = unique(points, 'rows'); % remove rows having same X-Y values
            
            isInside = arrayfun(@(x)CheckIsInside( avgFlow(t,:) , cropShot.*ratio , points(x,:) , maxSize ) , 1:size( points(:,1) ));
            outsiders = points(find(isInside == 0),:);
            extremes = [ min(outsiders(:,1)) max(outsiders(:,1)) ;
                         min(outsiders(:,2)) max(outsiders(:,2))];
                     
            % enlarge crop size for that shot if any point is left
            if ~(isempty(extremes))
                cropShot = ChangeCropSize( extremes , avgFlow(t,:) , cropShot , ratio , maxSize);
            end
        end
        crop = [crop; repmat(cropShot , shotEnd-shotStart+1 , 1)];
    end
    crop = [crop; crop(length(crop))];
end

function [crop] = ChangeCropSize( extremes , center , minSize , ratio , maxSize)
    X = extremes(1,:);
    Y = extremes(2,:);
    crop = minSize;
    
    if X(1,2) > center(1,1) + (crop*ratio(1,1))
        crop = ceil( (double(X(1,1)) - center(1,1)) / ratio(1,1)); 
    end
    if Y(1,2) > center(1,2) + (crop*ratio(1,2))
       crop = ceil( (double(Y(1,2)) - center(1,2)) / ratio(1,2)); 
    end
    if X(1,1) < center(1,1) - (crop*ratio(1,1))
       crop = floor( (center(1,1) - double(X(1,1))) / ratio(1,1)); 
    end
    if Y(1,2) < center(1,2) - (crop*ratio(1,2))
       crop = floor( (center(1,2) - double(Y(1,2))) / ratio(1,2)); 
    end 
    if (  center(1,1) - crop*ratio(1,1) < 0  ||...
            center(1,2) - crop*ratio(1,2) < 0 ||...
            crop*ratio(1,1) + center(1,1) > maxSize(1,1) ||...
            crop*ratio(1,2) + center(1,2) > maxSize(1,2) )
            
        crop = minSize;
    end
end

function [bool] = CheckIsInside(center , minSize , point , maxSize)
    bool = 1;
    if (point(1,1) < 0) || (point(1,2) < 0) ||...
            (point(1,1) > maxSize(1,1)) || (point(1,2) > maxSize(1,2))
        return;
    end
    if (( center(1,1) + minSize(1,1) < point(1,1) ) ||...
        ( center(1,1) - minSize(1,1) > point(1,1) ) ||...
        ( center(1,2) + minSize(1,2) < point(1,2) ) ||...
        ( center(1,2) - minSize(1,2) > point(1,2) ) )
        bool = 0;
    end
end








