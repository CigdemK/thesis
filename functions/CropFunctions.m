function F = CropFunctions
    F.GetWindowSize = @GetWindowSize;
    F.CreateWindow = @CreateWindow;
    F.RemoveBlackBars = @RemoveBlackBars;
end

% PUBLIC FUNCTIONS

function cropAll = GetWindowSize(avgFlow, impPts, shotBoundaries, newSize, maxSize)
 
    nrShots = length(shotBoundaries)-1;
    Util = UtilFunctions;
    
    cropAll = [];
    for k = 1:nrShots
        
        cropShot = zeros(shotBoundaries(k+1)-shotBoundaries(k),2);
        for t = shotBoundaries(k) : shotBoundaries(k+1)-1     

            % get all gaze pts in each frame
            indices = find( impPts( : , 2 ) == t );
            nrIndices = size(indices,1);
            
            if nrIndices == 0; continue; end;
                
            points = impPts(indices, 3:4);
            points = unique(points, 'rows'); % remove rows having same X-Y values
            
            isInside = Util.CheckIsInside(avgFlow(t,:), points, newSize);
            outsiders = points(isInside == 0,:);
            extremes = [ min(outsiders(:,1)) max(outsiders(:,1)) ;
                         min(outsiders(:,2)) max(outsiders(:,2))];
                     
            % enlarge crop size for that shot if any point is left
            if ~(isempty(extremes))
                cropShot(t,:) = ChangeCropSize( extremes , avgFlow(t,:) , newSize, maxSize);
            else
                cropShot(t,:) = newSize;
            end
            
        end
        
        cropShot = round(max(cropShot));
        cropAll = [cropAll; cropShot];
        
    end
    
end

function [croppedX,croppedY] = CreateWindow(cropSize, shotBoundaries, avg, maxSize)
    
    croppedX = [];
    croppedY = [];
    nrShots = length(shotBoundaries)-1;
    
    for i = 1:nrShots
        
        for k = shotBoundaries(i) : shotBoundaries(i+1)-1     
    
            minX = floor( avg( k , 1) - cropSize(i,1) /2) ;
            maxX = floor( avg( k , 1) + cropSize(i,1) /2) ;
            minY = floor( avg( k , 2) - cropSize(i,2) /2) ;
            maxY = floor( avg( k , 2) + cropSize(i,2) /2) ;

            if minX <= 0 
                minX = 1;
                maxX = cropSize(i,1) + 1 ;
            end
            if maxX > maxSize(1) 
                maxX = maxSize(1);
                minX = maxSize(1) - cropSize(i,1);
            end
            if minY <= 0 
                minY = 1;
                maxY = cropSize(i,2) + 1;
            end
            if maxY > maxSize(2) 
                maxY = maxSize(2) ;
                minY = maxSize(2) - cropSize(i,2);
            end
            croppedX = [croppedX ; floor([minX maxX]) ];
            croppedY = [croppedY ; floor([minY maxY]) ];
        end
    end
    croppedX = [croppedX ; croppedX(end,:) ];
    croppedY = [croppedY ; croppedY(end,:) ];
    
end

function [frames] = RemoveBlackBars(frames)

    [vHeight, vWidth, vDepth, nFrames] = size(frames);

    firstFrame = frames(:,:,:,1);
    
    % Detect the rows to be removed
    indices = [];
    for j = 1:vHeight
        if (length(firstFrame(firstFrame(j,:,:)<20)) == vWidth*vDepth)
            indices = [indices;j];   
        end
    end
    
    % Remove rows from all frames
    frames(indices,:,:,:) = [];
    
    indices = [];
    % Detect columns to be removed
    for j = 1:vWidth
        if (length(firstFrame(firstFrame(:,j,:)<20)) == vHeight*vDepth)
            indices = [indices;j];   
        end
    end
    
    % Remove columns from all frames
    frames(:,indices,:,:) = [];

end

% PRIVATE FUNCTIONS

function cropSize = ChangeCropSize( extremes , center , newSize, maxSize)

    X = extremes(1,:);
    Y = extremes(2,:);
    cropRatio = newSize(1)/newSize(2); % must preserve ratio

    XDiff = abs(X - center(2)) - newSize(1)/2;
    YDiff = abs(Y - center(1)) - newSize(2)/2;   
    
    XDiff = XDiff(XDiff>0);
    YDiff = YDiff(YDiff>0);
    
    cropSize(1) = newSize(1)/2 + max([XDiff 0]);
    cropSize(2) = newSize(2)/2 + max([YDiff 0]);
    
    loopcounter = 0;
    while(true)
        
        % Equalize the crop ratio
        if cropSize(1)/cropSize(2) > cropRatio
            cropSize(1) = cropSize(2)*cropRatio;
        else
            cropSize(2) = cropSize(1)/cropRatio;
        end

        % check boundary violations
        tooLargeX = center(2) + cropSize(1) > maxSize(1);
        tooLargeY = center(1) + cropSize(2) > maxSize(2);
        tooSmallX = center(2) - cropSize(1) < 1;
        tooSmallY = center(1) - cropSize(2) < 1;

        if tooSmallX;  cropSize(1) = center(2); end;
        if tooSmallY;  cropSize(2) = center(1); end;
        if tooLargeX;  cropSize(1) = maxSize(1) - center(2); end;
        if tooLargeY;  cropSize(2) = maxSize(2) - center(1); end;
        
        if abs(cropSize(1)/cropSize(2) - cropRatio)<0.00001; break; end;
        if loopcounter > 100
            break;
        else
            loopcounter = loopcounter+1;
        end
    end

    cropSize = cropSize .*2;
    
%     if X(1,2) > center(1,1) + minSize(1,1)
%         cropSize = ceil( double(X(1,2)) - center(1,1)); 
%     end
%     if Y(1,2) > center(1,2) + minSize(1,2)
%        cropSize = ceil( double(Y(1,2)) - center(1,2)); 
%     end
%     if X(1,1) < center(1,1) - minSize(1,1)
%        cropSize = floor( center(1,1) - double(X(1,1))); 
%     end
%     if Y(1,2) < center(1,2) - minSize(1,2)
%        cropSize = floor( center(1,2) - double(Y(1,2))); 
%     end 
%     if (  minSize(1,1) + center(1,1) > maxSize(1,1) ||...
%           minSize(1,2) + center(1,2) > maxSize(1,2) )
%         cropSize = min((maxSize(1,2)-center(1,2)), (maxSize(1,1)-center(1,1)));
%     end
%     if ( center(1,1) - minSize(1,1) < 0  ||...
%          center(1,2) - minSize(1,2) < 0 )
%         cropSize = max(center(1,2),center(1,1));
%     end
end

