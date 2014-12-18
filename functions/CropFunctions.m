function F = CropFunctions
    F.EstimateCropWindowSize = @EstimateCropWindowSize;
    F.CreateCropWindow = @CreateCropWindow;
    F.RemoveBlackBars = @RemoveBlackBars;
end

% PUBLIC FUNCTIONS

function [crop] = EstimateCropWindowSize(avgFlow, eyes , shotBoundaries ,minSize,maxSize)
 
    crop = [];
    nrShots = length(shotBoundaries)-1;
    Util = UtilFunctions;

    for k = 1:nrShots
        shotStart = shotBoundaries(k);
        shotEnd = shotBoundaries(k+1)-1;
        shotFlow = avgFlow( shotStart:shotEnd , :);
        
        for t = shotStart : shotEnd     
            
            % get all gaze pts in each frame
            indices = find( eyes( : , 2 ) == t );
            nrIndices = length(indices);
            Fx = eyes( indices( 1:nrIndices ) , 4 );
            Fy = eyes( indices( 1:nrIndices ) , 3 );
            points = [Fx Fy];
            points = unique(points, 'rows'); % remove rows having same X-Y values
            
            isInside = arrayfun(@(x)Util.CheckIsInside( avgFlow(t,:) , minSize , points(x,:) , maxSize ) , 1:size( points(:,1) ));
            outsiders = points(find(isInside == 0),:);
            extremes = [ min(outsiders(:,1)) max(outsiders(:,1)) ;
                         min(outsiders(:,2)) max(outsiders(:,2))];
                     
            % enlarge crop size for that shot if any point is left
            if ~(isempty(extremes))
                cropShot(t) = ChangeCropSize( extremes , avgFlow(t,:) , minSize , maxSize);
            end
        end
        cropShot = max(cropShot);
        crop = [crop; repmat(cropShot , shotEnd-shotStart+1 , 1)];
    end
    crop = [crop; crop(length(crop),:)];
end

function [croppedX,croppedY] = CreateCropWindow( limX , limY , avg , vidWidth , vidHeight)
    
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
        croppedX = [croppedX ; round([minX maxX]) ];
        croppedY = [croppedY ; round([minY maxY]) ];
    end
    
end

function [frames] = RemoveBlackBars(frames)

    [vHeight, vWidth, vDepth, nFrames] = size(frames);
    indices = [];
    
    % Detect the rows to be removed
    for j = 1:vHeight
        line = frames(:,:,:,1);
        if (length(line(line(1,:,:)<20)) == vWidth*vDepth)
            indices = [indices;j];   
        end
    end
    
    % Remove rows from all frames
    for n = 1:nFrames
        frames(indices,:,:,n) = [];
    end
end

% PRIVATE FUNCTIONS

function [crop] = ChangeCropSize( extremes , center , minSize , maxSize)
    X = extremes(1,:);
    Y = extremes(2,:);
    crop = minSize;
    
    if X(1,2) > center(1,1) + minSize(1,1)
        crop = ceil( double(X(1,2)) - center(1,1)); 
    end
    if Y(1,2) > center(1,2) + minSize(1,2)
       crop = ceil( double(Y(1,2)) - center(1,2)); 
    end
    if X(1,1) < center(1,1) - minSize(1,1)
       crop = floor( center(1,1) - double(X(1,1))); 
    end
    if Y(1,2) < center(1,2) - minSize(1,2)
       crop = floor( center(1,2) - double(Y(1,2))); 
    end 
    if (  minSize(1,1) + center(1,1) > maxSize(1,1) ||...
          minSize(1,2) + center(1,2) > maxSize(1,2) )
        crop = min((maxSize(1,2)-center(1,2)), (maxSize(1,1)-center(1,1)));
    end
    if ( center(1,1) - minSize(1,1) < 0  ||...
         center(1,2) - minSize(1,2) < 0 )
        crop = max(center(1,2),center(1,1));
    end
end