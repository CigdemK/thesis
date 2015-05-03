function F = WANGFunctions
    F.ScaleStretch = @ScaleStretch;
    F.PerFrameOptimization = @PerFrameOptimization;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Public Functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [retargetedFrames,gridX,gridY] = PerFrameOptimization(frames, newsize,videoName)

    warning ('off','all');
    
%   Remove black bars and adjust the size parameters
    Cropper = CropFunctions;
    frames = Cropper.RemoveBlackBars(frames);
    [vidHeight,vidWidth,~,nFrames] = size(frames);
    
% Step 1: Apply Scale&Stretch to all frames
    meshSize = 20;
    originalFrames = zeros(vidHeight,vidWidth,3,nFrames);
    scaledFrames = zeros(newsize(2),newsize(1),3,nFrames);
    gridX = zeros(ceil(vidHeight/meshSize),ceil(vidWidth/meshSize),nFrames);
    gridY = zeros(ceil(vidHeight/meshSize),ceil(vidWidth/meshSize),nFrames);

    disp('Scaling for the first time');
    for k = 1:nFrames
        originalFrames(:,:,:,k) = frames(:,:,:,k);
%         [tmp, gridX(:,:,k), gridY(:,:,k)]= ScaleStretch(frames(:,:,:,k), ...
%             [newsize(2) newsize(1)],[meshSize meshSize]);
%         scaledFrames(:,:,:,k) = im2double(tmp);
    end
    
% save(['/home/gul.varol/cigdem/AVIClips/scaled/' videoName '_scaledFrames.mat'],'scaledFrames','gridX','gridY');
% Step 2: Calculate pathlines and adjacencies
%     tmp=regexp(videoName,'\','split');
% 
%     load([tmp{end} '_Frames.mat_scaledFrames.mat'],'scaledFrames','gridX','gridY');
%     disp('Compute pathlines');
% %     originalPathlines   = GetPathlines( originalFrames, [vidHeight/20 vidWidth/20]);
%     load('originalPathlines.mat','originalPathlines')
%     scaledPathlines     = GetWarpedPathlines(originalPathlines,gridX,gridY,[vidHeight,vidWidth]);
%     nPathlines = size(originalPathlines,2);
%     
%     originalAdjacencies = CalculateAdjacencies(originalPathlines); 
%     save([tmp{end} '_pathlines.mat'] , 'originalPathlines','scaledPathlines','originalAdjacencies');  nPathlines = size(originalPathlines,2); 
  
%     load(['/home/gul.varol/cigdem/AVIClips/pathlines/' videoName '_pathlines.mat'] , 'originalPathlines','scaledPathlines','originalAdjacencies');  nPathlines = size(originalPathlines,2); 
%     
%     S = ones(2,nPathlines) ;%+ 0.1* rand(2,nCorrespondances);
%     t = zeros(2,nPathlines) ;%+ 0.1* rand(2,nCorrespondances);
%     Sij = ones(4,nPathlines) ;%+ 0.1* rand(4,nCorrespondances);
% 
%     pstart = [originalPathlines.start];
%     plength = [originalPathlines.pend] - pstart+1;

 %  Minimize the error between the corresponding pathlines 
%     disp('Entering S loop');
%     for k = 1:3
%         for i = 1:nPathlines
%             adjacencyIndices =  originalAdjacencies(originalAdjacencies(:,1) == i,2);
%             nAdjacencies = size(adjacencyIndices,1);
%             
%             if nAdjacencies == 0; continue; end;
%             
%             % Prepare inputs for error function
%             adjacencyFrames = [];
%             adjacencies = []; % [x1 y1 1 x2 y2 1 (x1-x2) (y1-y2)]'
%             for j = 1:nAdjacencies
%                 
%                 % Find which pathline starts first
%                 start1 = originalPathlines(i).start;
%                 start2 = originalPathlines(adjacencyIndices(j)).start;
%                 end1 = originalPathlines(i).pend;
%                 end2 = originalPathlines(adjacencyIndices(j)).pend;
%                 currentStart = max(start1,start2);
%                 currentEnd = min(end1,end2);
%                 currentLength = currentEnd - currentStart + 1;
%                 ind1 = currentStart - start1 + 1: currentEnd - start1 + 1;
%                 ind2 = currentStart - start2 + 1: currentEnd - start2 + 1;
%                       
%                 adjacencyFrames = [adjacencyFrames,repmat(j,[1 currentLength])];
%                 adjacencies = [adjacencies,...
%                     [reshape(originalPathlines(i).x(ind1),[1 currentLength]);...
%                     reshape(originalPathlines(i).y(ind1),[1 currentLength]);...
%                     ones(1,currentLength);...
%                     reshape(originalPathlines(adjacencyIndices(j)).x(ind2),[1 currentLength]);...
%                     reshape(originalPathlines(adjacencyIndices(j)).y(ind2),[1 currentLength]);...
%                     ones(1,currentLength);...
%                     reshape(originalPathlines(i).x(ind1) - originalPathlines(adjacencyIndices(j)).x(ind2),...
%                         [1 currentLength]);...
%                     reshape(originalPathlines(i).y(ind1) - originalPathlines(adjacencyIndices(j)).y(ind2),...
%                         [1 currentLength])]];
%             end
% 
%             % Omege function calculates and outputs the error
%              f = @(S)OmegaFunctions(S, [reshape(originalPathlines(i).x,[1 numel(originalPathlines(i).x)]);...
%                                         reshape(originalPathlines(i).y,[1 numel(originalPathlines(i).x)]);...
%                                         ones(1,numel(originalPathlines(i).x))],...
%                                         [reshape(scaledPathlines(i).x,[1 numel(originalPathlines(i).x)]);...
%                                         reshape(scaledPathlines(i).y,[1 numel(originalPathlines(i).x)])],...
%                                         adjacencies, adjacencyFrames ) ;
% 
%             % Prepare inputs for error function
%             currentS = [S(:,i);t(:,i)]; % [Six Siy tix tiy ...]'
%             for j = 1:nAdjacencies
%                 currentS = [currentS;... %[... Sjx Sjy tjx tjy Sijxx Sijyx Sijxy Sijyy]'
%                    (-1)*S(:,adjacencyIndices(j));...
%                    (-1)*t(:,adjacencyIndices(j));...
%                    (-1)*Sij(:,adjacencyIndices(j))];
%             end
% 
%             % Run the error function optimization
%             option = optimset('Display','off','LargeScale','off','Algorithm','interior-point');
%             optimalS = fminunc(f,currentS,option);
%             
%             S(:,i) = optimalS(1:2);
%             t(:,i) = optimalS(3:4);
%             for j = 1:nAdjacencies
%                 currInd = 5+8*(j-1);
%                 Sij(:,adjacencyIndices(j)) = optimalS(currInd+4:currInd+7);
%             end  
%         end
%         disp('one gone in the fifth');
%         save([videoName '_S.mat'],'S','t','Sij');
%     end
    
%     save([videoName '_S.mat'],'S');
%     load(['/home/gul.varol/cigdem/AVIClips/s/' videoName '_S.mat']);
% Step 3: Optimize Scaling with the Pathlines  
load(['/home/gul.varol/cigdem/AVIClips/pathlines/' videoName '_pathlines.mat'] , 'originalPathlines','scaledPathlines','originalAdjacencies');  nPathlines = size(originalPathlines,2); 
load(['/home/gul.varol/cigdem/AVIClips/scaled/' videoName '_scaledFrames.mat'],'scaledFrames','gridX','gridY');
load(['/home/gul.varol/cigdem/AVIClips/s/' videoName '_S.mat']);
%     Calculate positional constraints

    pathlineLocations = struct('x',[],'y',[]);
    for i = 1:nPathlines
        pathlineLocations(i).x = originalPathlines(i).x .* S(1,i);
        pathlineLocations(i).y = originalPathlines(i).y .* S(2,i);
    end
    
    mystart = [originalPathlines.start];
    myend = [originalPathlines.pend];
    
    disp('Start rescaling')
    for k = 1:nFrames
        
        currentPathlines = mystart<=k & myend>=k; 
        currentIndices = k - mystart + 1;
        
        passingPathlines = [];
        for i = 1:nPathlines
            if currentPathlines(i)
                passingPathlines = [passingPathlines, ...
                    [pathlineLocations(i).x(currentIndices(i)); ...
                    pathlineLocations(i).y(currentIndices(i))]];
            end
        end
        
        passingPathlines(1,passingPathlines(1,:)>vidWidth) = vidWidth;
        passingPathlines(2,passingPathlines(2,:)>vidHeight) = vidHeight;
        passingPathlines(passingPathlines<1) = 1;
        
        if ~isempty(passingPathlines)
            [tmp, ~, ~]= ScaleStretch(frames(:,:,:,k),[newsize(2) newsize(1)],[20 20],passingPathlines);
            scaledFrames(:,:,:,k) = im2double(tmp);
        end
    end
    
    retargetedFrames = scaledFrames;
end

function [imgOut, meshXNew , meshYNew]= ScaleStretch(imgIn, newsize, meshsize, constraints)

    [height, width, ~] = size(imgIn);

    % Compute significance map
    Wb = simpsal(imgIn);
    [Gx,Gy] = imGradient(double(imgIn(:,:,1)));
    Wa = sqrt(Gx.^2 + Gy.^2);
    W = mat2gray(Wa.*imresize(Wb,[height, width]));

    % Construct the mesh grid
    [meshX,meshY] = meshgrid(1:meshsize(2):width,1:meshsize(1):height);
    meshX(:,end) = width;
    meshY(end,:) = height;
        
    Y = size(meshY,1);
    X = size(meshX,2);
    numberOfQuads = ( Y - 1 ) * ( X - 1 );  

    % Calculate grid weights
    Wq = zeros(numberOfQuads*4,1);
    for w = 1:X-1
        for h = 1:Y-1
            y = meshY(h,w):meshY(h+1,w);
            x = meshX(h,w):meshX(h,w+1);
            Wq(4*((w-1)*(Y-1)+h)-3) =  mean( mean( W( y , x ) ) );
            Wq(4*((w-1)*(Y-1)+h)-2) =  mean( mean( W( y , x ) ) );
            Wq(4*((w-1)*(Y-1)+h)-1) =  mean( mean( W( y , x ) ) );
            Wq(4*((w-1)*(Y-1)+h))   =  mean( mean( W( y , x ) ) );
        end 
    end
    Wq = Wq / (max(max(Wq)));

    % New meshes are the startig points of the iteration
    meshXNew = ceil(meshX/(meshX(end:end)/newsize(2)));
    meshYNew = ceil(meshY/(meshY(end:end)/newsize(1)));
    meshXNew(:,end) = newsize(2);
    meshYNew(end,:) = newsize(1);
    meshXOld = zeros(Y,X);
    meshYOld = zeros(Y,X);
    
    % Create a mask to convert vertices to vertex differences. Vertex
    % differences are specific to each quad, so there are in total
    % numberOfQuads*4 vertex differences (correnpond to edges). These
    % matrices are precomputed to be used in all iterations.
    maskDiff = VertexDiffMask(X,Y);
    originalVerticesX = maskDiff*meshX(:);
    originalVerticesY = maskDiff*meshY(:);
    
    if exist('constraints','var') == 1
        maskVert = VertexMask(X,Y);
        closestVerticesX = maskVert*constraints(1,:)';
        closestVerticesY = maskVert*constraints(2,:)';
    end
    % Set limits of fmincon so that the output image will be
    % rectangular
    options = optimset('Algorithm','interior-point','Display','off');
    AeqX = zeros(Y*X);
    beqX = zeros(Y*X,1);
    AeqY = zeros(Y*X);
    beqY = zeros(Y*X,1);
    lbX = ones(Y*X,1);
    ubX = ones(Y*X,1)*newsize(2);
    ubY = ones(Y*X,1)*newsize(1);
    
    for i = 0:Y-1 
        AeqX(end-i,end-i) = 1;
        AeqX(i+1,i+1) = 1;
        beqX(end-i) = newsize(2);
        beqX(i+1) = 1;
    end
    for i = 0:X-1 
        AeqY(i*Y+1,i*Y+1) = 1;
        AeqY(i*Y+Y,i*Y+Y) = 1;
        beqY(i*Y+Y) = newsize(1);
        beqY(i*Y+1) = 1;
    end

    oldErr = 10000000;
    t = 1;
    while (1)

        [res,err] = CheckIfDone(meshXOld , meshYOld  , meshXNew , meshYNew );
        if res == 1
            break;
        end
        if oldErr < err
            meshXOld = meshXOld * 0.3 + meshXNew*0.7;
            meshYOld = meshYOld * 0.3 + meshYNew*0.7;
            t = t + 1;
        else
            oldErr = err;
            meshXOld = meshXNew;
            meshYOld = meshYNew;
        end
        
        if t == 5
            break;
        end
        
        % Calculate S and l with the current estimates of V'
        S = zeros(numberOfQuads*4,1);
        l = zeros(numberOfQuads*4,1);
        for k = 1:numberOfQuads
            
            yInd = rem(k-1,Y-1)+1;
            xInd = floor((k-1)/(Y-1))+1;
            currentQuad = GetQuad_2Dim(meshX, meshY, xInd, yInd);
            currentQuadNew = GetQuad_2Dim(meshXOld, meshYOld, xInd, yInd);
    
            top = 0;
            bottom = 0;
            for i = 1:4 %for each edge of a quad
                
                ind = mod(i-1,4)+1;
                nextInd = mod(i,4)+1;
                
                top = top +  (currentQuad(:,ind) - currentQuad(:,nextInd))' * ...
                    (currentQuadNew(:,ind) - currentQuadNew(:,nextInd));
                bottom = bottom +  (currentQuad(1,ind) - currentQuad(1,nextInd))^2 + ...
                    (currentQuad(2,ind) - currentQuad(2,nextInd))^2;
                
                l((k-1)*4+i) = ( (currentQuadNew(1,ind) - currentQuadNew(1,nextInd))^2 + ...
                    (currentQuadNew(2,ind) - currentQuadNew(2,nextInd))^2 ) /...
                    ( (currentQuad(1,ind) - currentQuad(1,nextInd))^2 + ...
                    (currentQuad(2,ind) - currentQuad(2,nextInd))^2 );
                
            end
            S(k*4-3) = top / bottom;
            S(k*4-2) = top / bottom;
            S(k*4-1) = top / bottom;
            S(k*4)   = top / bottom;
  
        end

        if exist('constraints','var') == 1
            
            FX = @(Vnew)sum(((maskDiff*Vnew).^2) .* (Wq + 1) + ...
                             (originalVerticesX.^2) .* (Wq .* (S.^2) + (l.^2)) + ...
                             (maskDiff*Vnew).*(originalVerticesX) .* ((Wq .* S + l) * (-2)) + ...
                             abs(maskDiff*Vnew .* Wq - closestVerticesX));
            FY = @(Vnew)sum(((maskDiff*Vnew).^2) .* (Wq + 1) + ...
                             (originalVerticesY.^2) .* (Wq .* (S.^2) + (l.^2)) + ...
                             (maskDiff*Vnew).*(originalVerticesY) .* ((Wq .* S + l) * (-2)) + ...
                             abs(maskDiff*Vnew .* Wq - closestVerticesY));
        else
            FX = @(Vnew)sum(((maskDiff*Vnew).^2) .* (Wq + 1) + ...
                             (originalVerticesX.^2) .* (Wq .* (S.^2) + (l.^2)) + ...
                             (maskDiff*Vnew).*(originalVerticesX) .* ((Wq .* S + l) * (-2)));
            FY = @(Vnew)sum(((maskDiff*Vnew).^2) .* (Wq + 1) + ...
                             (originalVerticesY.^2) .* (Wq .* (S.^2) + (l.^2)) + ...
                             (maskDiff*Vnew).*(originalVerticesY) .* ((Wq .* S + l) * (-2)));
        end
        
        meshXNew = fmincon( FX ,meshXOld(:) ,[] ,[] ,AeqX ,beqX ,lbX ,ubX ,[], options );
        meshYNew = fmincon( FY ,meshYOld(:) ,[] ,[] ,AeqY ,beqY ,lbX ,ubY ,[], options );
        
        meshXNew = reshape(meshXNew, [Y X] );
        meshYNew = reshape(meshYNew, [Y X] );

    end

    % Warp image with the new mesh
    meshXNew = reshape(meshXOld,Y,X);
    meshYNew = reshape(meshYOld,Y,X);

    imgOut = WarpImageWithMesh(imgIn , meshXNew , meshYNew);    
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Private Functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Functions for PerFrameOptimization

function [pathlines] = GetPathlines(frames , meshSize, opticalFlow)
    
    % Calculate optical flow
    if nargin < 3
        Motion = MotionFunctions;
        [Vx, Vy] = Motion.MyOpticalFlow(frames);
    else
        Vx = opticalFlow(:,:,:,1);
        Vy = opticalFlow(:,:,:,2);
    end

    [vidHeight , vidWidth, ~, nFrames] = size(frames)
    [vidHeightO , vidWidthO, nFramesO] = size(Vx)
    
    % Initialize mesh grids
    [ initMeshX , initMeshY ] = meshgrid( 1:floor(vidWidth/meshSize(2)):vidWidth, ...
                                          1:floor(vidHeight/meshSize(1)):vidHeight );
                                      
    X = size(initMeshX,2);
    Y = size(initMeshX,1);
    
    meshX = zeros(Y,X,nFrames);
    meshY = zeros(Y,X,nFrames);
    
    meshX(:,:,1) = initMeshX;
    meshY(:,:,1) = initMeshY;

    % Preallocation
    maxNrOfPathlines = Y*X*nFrames;
    [ xReg , yReg ] = meshgrid(1:vidWidth,1:vidHeight);
    pathlinesStart = ones(Y,X);
    pathlines(maxNrOfPathlines) = struct('start',[],'pend',[],'x',[],'y',[],'xind',[],'yind',[]);
    
    j = 1;
    for k = 2:nFrames
               
        % Find optical flow values on mesh scale
        xPts = interp2(xReg , yReg , Vx(:,:,k-1) , meshX(:,:,k-1) , meshY(:,:,k-1) );
        yPts = interp2(xReg , yReg , Vy(:,:,k-1) , meshX(:,:,k-1) , meshY(:,:,k-1) );
        
        % Add the current optical flow to the prev. vertex pts
        meshX(:,:,k) = meshX(:,:,k-1) + xPts;
        meshY(:,:,k) = meshY(:,:,k-1) + yPts;
        
        % Find terminated paths that exceeds frame size
        [indXr, indXc] = find( meshX(:,:,k)<1 | meshX(:,:,k)>vidWidth );
        [indYr, indYc] = find( meshY(:,:,k)<1 | meshY(:,:,k)>vidHeight );
        
        indices = unique([[indXr; indYr], [indXc; indYc]],'rows');        

        for i = 1:size(indices,1)

            % Add pathline terminated paths to the path list
            pathlines(j).start = pathlinesStart(indices(i,1),indices(i,2));
            pathlinesStart(indices(i,1),indices(i,2)) = k;
            pathlines(j).x = meshX(indices(i,1),indices(i,2),pathlines(j).start:k-1);
            pathlines(j).y = meshY(indices(i,1),indices(i,2),pathlines(j).start:k-1);
            pathlines(j).x = round(pathlines(j).x(:));
            pathlines(j).y = round(pathlines(j).y(:));
            pathlines(j).xind = indices(i,1);
            pathlines(j).yind = indices(i,2);
            pathlines(j).pend = pathlines(j).start + size(pathlines(j).x,1)-1;
            j = j + 1;

            % Add new seeds
            meshX(indices(i,1),indices(i,2),k) = initMeshX(indices(i,1),indices(i,2));
            meshY(indices(i,1),indices(i,2),k) = initMeshY(indices(i,1),indices(i,2));

        end
        
    end
    
    % Add non terminated pathlines at the end

    for i = 1:Y
        for k = 1:X
            pathlines(j).start = pathlinesStart(i,k);
            pathlines(j).x = meshX(i,k,pathlines(j).start:end);
            pathlines(j).y = meshY(i,k,pathlines(j).start:end);
            pathlines(j).x = round(pathlines(j).x(:));
            pathlines(j).y = round(pathlines(j).y(:));
            pathlines(j).xind = k;
            pathlines(j).yind = i;
            pathlines(j).pend = pathlines(j).start + size(pathlines(j).x,1)-1;
            j = j + 1;
        end
    end
    
    % remove all zero lines
    pathlines(j:end) = [];
    
end

function [adjacencies] = CalculateAdjacencies(pathlines)

    adjacencies = [];
    nPathlines = size(pathlines,2);
    
    start = [pathlines.start];
    pend = [pathlines.pend];
    xind = [pathlines.xind];
    yind = [pathlines.yind];
    
    for i = 1:nPathlines

        inds = start<=start(i) & pend>=pend(i);
        xDiff = xind.*inds - xind(i);
        yDiff = yind.*inds - yind(i);
        
        candidates = find(sum([abs(xDiff'),abs(yDiff')],2)==1);
        adjacencies = [adjacencies; [ repmat(i,[numel(candidates),1]) candidates ] ] ;
       
    end
    adjacencies = unique(sort(adjacencies,2),'rows');
    
end

function scaledPathlines = GetWarpedPathlines(originalPathlines,gridX,gridY,originalSize)

    nPathlines = size(originalPathlines,2);
    [Y,X,nFrames] = size(gridX);
    scaledSizeX = max(gridX(:));
    scaledSizeY = max(gridY(:));
    
    scaledPathlines(nPathlines) = struct('start',[],'pend',[],'x',[],'y',[]);
    
%     [xBig, yBig, zBig] = meshgrid(1:originalSize(2),1:originalSize(1),1:nFrames);
    
%     meshSizeX = 1:ceil(scaledSizeX/X):round(scaledSizeX);
%     meshSizeY =  1:ceil(scaledSizeY/Y):round(scaledSizeY);
%     if numel(meshSizeX)<X
%         meshSizeX(end+1) = X;
%     end
%     if numel(meshSizeY)<Y
%         meshSizeY(end+1) = Y;
%     end
    
%     [xSmall, ySmall, zSmall] = meshgrid(meshSizeX, meshSizeY,1:nFrames);
%     xSmall(:,end,:) = scaledSizeX;
%     ySmall(end,:,:) = scaledSizeY;
    
    meshSizeX = 1:scaledSizeX/originalSize(2):scaledSizeX;
    meshSizeY = 1:scaledSizeY/originalSize(1):scaledSizeY;
    if numel(meshSizeX)<originalSize(2)
        meshSizeX(end+1) = scaledSizeX;
    end
    if numel(meshSizeY)<originalSize(1)
        meshSizeY(end+1) = scaledSizeY;
    end
    [xBig, yBig, zBig] = meshgrid(meshSizeX,meshSizeY,1:nFrames);
    
    meshSizeX = 1:ceil(originalSize(2)/X):round(originalSize(2));
    meshSizeY = 1:ceil(originalSize(1)/Y):round(originalSize(1));
    if numel(meshSizeX)<X
        meshSizeX(end+1) = originalSize(2);
    end
    if numel(meshSizeY)<Y
        meshSizeY(end+1) = originalSize(1);
    end
    
    [xSmall, ySmall, zSmall] = meshgrid(meshSizeX, meshSizeY,1:nFrames);
    xSmall(:,end,:) = originalSize(2);
    ySmall(end,:,:) = originalSize(1);

    newCoordsX = interp3( xSmall, ySmall, zSmall, gridX, xBig, yBig, zBig, 'spline');
    newCoordsY = interp3( xSmall, ySmall, zSmall, gridY, xBig, yBig, zBig, 'spline');
    
    start = [originalPathlines.start];
    plength = [originalPathlines.pend] - start + 1;
    [scaledPathlines.start] = originalPathlines.start;
    [scaledPathlines.pend] = originalPathlines.pend;
    
    for k = 1:nPathlines
             
        theX = zeros(plength(k),1);
        theY = zeros(plength(k),1);
        for t = 1:plength(k)
            theX(t) = newCoordsX(originalPathlines(k).y(t),originalPathlines(k).x(t),start(k)+t-1);
            theY(t) = newCoordsY(originalPathlines(k).y(t),originalPathlines(k).x(t),start(k)+t-1);
        end
        
        scaledPathlines(k).x = theX;
        scaledPathlines(k).y = theY;
   
    end
    
end

function y = OmegaFunctions(allInputs, originalPath, scaledPath, adjacencies, adjacencyFrames)
% allInputs = [Six Siy tix tiy, 
             % Sjx Sjy tjx tjy Sijxx Sijxy Sijyy Sijyx]'; 
             % (repeat the second line for all adjacencies)
             
    [~, nFrames] = size(originalPath);
    nAdjacencies = ( size(allInputs,1) - 4 ) / 8;

    SD = [allInputs(1) 0 allInputs(3); 0 allInputs(2) allInputs(4)];
    for i = 1: nAdjacencies
        currInd = 5+8*(i-1);
        SP(:,:,i) = [allInputs(currInd) 0   allInputs(currInd+2) allInputs(currInd+4) allInputs(currInd+5);
                     0 allInputs(currInd+1) allInputs(currInd+3) allInputs(currInd+6) allInputs(currInd+7)];
    end
    
    SP = [repmat(SD,[1,1,nAdjacencies]) , SP];

    omegaD = 0;
    for i = 1:nFrames
        currentError =  SD * originalPath(:,i) - scaledPath(:,i);
        omegaD = omegaD + sqrt((currentError(1)^2)+(currentError(2)^2));
    end

    omegaP = 0;
    for k = 1: size(adjacencies,2)
        currentAdj = adjacencyFrames(k);
        currentError =  SP(:,:,currentAdj) * adjacencies(:,k);
        omegaP = omegaP + sqrt((currentError(1)^2)+(currentError(2)^2));
    end 
    
    y = omegaD * 0.5 + omegaP;
end


% % Functions for ScaleStretch

function binaryMatrix =  VertexDiffMask(X,Y)

    numberOfVertices = X*Y;
    numberOfQuads = (X-1)*(Y-1);
    
    vector = zeros(4,numberOfVertices);
    vector(1,1) = 1;     vector(2,Y+1) =  1;   vector(3,Y+2) = 1;    vector(4,2) = 1;
    vector(1,Y+1) = -1;  vector(2,Y+2) = -1;   vector(3,2)   = -1;   vector(4,1) = -1;
    
    binaryMatrix = zeros(numberOfQuads*4,numberOfVertices);
    for x = 1:X-1
        for y = 1:Y-1
            binaryMatrix(4*((x-1)*(Y-1)+y)-3,:) = vector(1,:);
            binaryMatrix(4*((x-1)*(Y-1)+y)-2,:) = vector(2,:);
            binaryMatrix(4*((x-1)*(Y-1)+y)-1,:) = vector(3,:);
            binaryMatrix(4*((x-1)*(Y-1)+y)  ,:) = vector(4,:);

            vector(1,:) = circshift(vector(1,:),[-1 1]);
            vector(2,:) = circshift(vector(2,:),[-1 1]);
            vector(3,:) = circshift(vector(3,:),[-1 1]);
            vector(4,:) = circshift(vector(4,:),[-1 1]);
        end
        vector(1,:) = circshift(vector(1,:),[-1 1]);
        vector(2,:) = circshift(vector(2,:),[-1 1]);
        vector(3,:) = circshift(vector(3,:),[-1 1]);
        vector(4,:) = circshift(vector(4,:),[-1 1]);
    end
end

function binaryMatrix =  VertexMask(X,Y)

    numberOfVertices = X*Y;
    numberOfQuads = (X-1)*(Y-1);
    
    vector = zeros(4,numberOfVertices);
    vector(1,1) = 1;     vector(2,2) =  1;   vector(3,Y+1) = 1;    vector(4,Y+2) = 1;
    
    binaryMatrix = zeros(numberOfQuads*4,numberOfVertices);
    for x = 1:X-1
        for y = 1:Y-1
            binaryMatrix(4*((x-1)*(Y-1)+y)-3,:) = vector(1,:);
            binaryMatrix(4*((x-1)*(Y-1)+y)-2,:) = vector(2,:);
            binaryMatrix(4*((x-1)*(Y-1)+y)-1,:) = vector(3,:);
            binaryMatrix(4*((x-1)*(Y-1)+y)  ,:) = vector(4,:);

            vector(1,:) = circshift(vector(1,:),[-1 1]);
            vector(2,:) = circshift(vector(2,:),[-1 1]);
            vector(3,:) = circshift(vector(3,:),[-1 1]);
            vector(4,:) = circshift(vector(4,:),[-1 1]);
        end
        vector(1,:) = circshift(vector(1,:),[-1 1]);
        vector(2,:) = circshift(vector(2,:),[-1 1]);
        vector(3,:) = circshift(vector(3,:),[-1 1]);
        vector(4,:) = circshift(vector(4,:),[-1 1]);
    end
end

function [result,error] = CheckIfDone( meshXOld , meshYOld  , meshX , meshY )

    X = abs(meshXOld(:) - meshX(:));
    Y = abs(meshYOld(:) - meshY(:));
    
    error = sum(sqrt(X.^2 + Y.^2));
    result = true;
    
    if size(X(X>0.5)) + size(Y(Y>0.5)) > 0
        result = false;
    end
    
    if ~isempty([find(round(meshX)<1); find(round(meshY)<1)])
        result = true;
    end
    
    if ~isempty([find(round(meshX) > round(meshX(end))); find(round(meshY) > round(meshY(end)))])
        result =  true;
    end
    
end

function quad = GetQuad_2Dim(  meshX, meshY, xInd, yInd )
    quad = [[ meshY( yInd   , xInd )   ; meshX( yInd   , xInd   )],...
            [ meshY( yInd+1 , xInd )   ; meshX( yInd+1 , xInd   )],...
            [ meshY( yInd   , xInd+1 ) ; meshX( yInd   , xInd+1 )],...
            [ meshY( yInd+1 , xInd+1 ) ; meshX( yInd+1 , xInd+1 )]];
end

function warpedImage = WarpImageWithMesh(imgIn , meshX , meshY )

    [ Y , X ] = size(meshX);
    n = meshX(end:end);
    m = meshY(end:end);
    
    imgIn = imresize(imgIn,[round(m),round(n)]);

    [xRegular,yRegular] = meshgrid(1:round(n),1:round(m));
    meshSizeX = 1 : ceil( n/X) : n;
    meshSizeY = 1 : ceil( m/Y) : m ;
    
    if numel(meshSizeX)<X
        meshSizeX(end+1) = n;
    end
    if numel(meshSizeY)<Y
        meshSizeY(end+1) = m;
    end
    [meshXRegular,meshYRegular] = meshgrid( meshSizeX , meshSizeY );
    meshXRegular(:,end) = n;
    meshYRegular(end,:) = m;
    
    newCoordsX = interp2( meshXRegular, meshYRegular, meshX, xRegular, yRegular, 'spline');
    newCoordsY = interp2( meshXRegular, meshYRegular, meshY, xRegular, yRegular, 'spline');
    
    warped1 = interp2( xRegular, yRegular, im2double(imgIn(:,:,1)), newCoordsX, newCoordsY, 'spline');
    warped2 = interp2( xRegular, yRegular, im2double(imgIn(:,:,2)), newCoordsX, newCoordsY, 'spline');
    warped3 = interp2( xRegular, yRegular, im2double(imgIn(:,:,3)), newCoordsX, newCoordsY, 'spline');
    
    warpedImage = uint8(cat(3,warped1,warped2,warped3).*255);

end

function PlotPathlines(originalPathlines,scaledPathlines,stepSize,...
    correspondances,originalAdjacencies,scaledAdjacencies , originalSize , scaledSize)

    subplot(3,1,1);
    title('Some Correspoding Pathlines'); hold on;
    legend('Scaled Video(Resized)','Original Video'); hold on;
    axis([1 originalSize(2) 1 originalSize(1)]); hold on;
    for i =1:stepSize:size(correspondances,1)
        x=scaledPathlines(correspondances(i,2),1,:).*(originalSize(2)/scaledSize(2));
        y=scaledPathlines(correspondances(i,2),2,:).*(originalSize(1)/scaledSize(1));
        set(gca,'YDir','reverse'); 
        plot(x(:),y(:),'rx'); 
        hold on;
    end
    
    for i =1:stepSize:size(correspondances,1)
        x=originalPathlines(correspondances(i,1),1,:); 
        y=originalPathlines(correspondances(i,1),2,:); 
        set(gca,'YDir','reverse');
        plot(x(:),y(:),'bx'); 
        hold on;
    end
    hold off;
    
    subplot(3,1,2);
    title('Some Adjacent Pathlines From Original Video'); hold on;
    axis([1 originalSize(2) 1 originalSize(1)]); hold on;
    for i =1:stepSize:size(originalAdjacencies,1)
        x=originalPathlines(originalAdjacencies(i,2),1,:); 
        y=originalPathlines(originalAdjacencies(i,2),2,:); 
        set(gca,'YDir','reverse');
        plot(x(:),y(:),'bx'); 
        hold on;
    end
    for i =1:stepSize:size(originalAdjacencies,1)
        x=originalPathlines(originalAdjacencies(i,1),1,:); 
        y=originalPathlines(originalAdjacencies(i,1),2,:); 
        plot(x(:),y(:),'rx'); 
        hold on;
    end
    hold off;
    
    subplot(3,1,3);
    title('Some Adjacent Pathlines From Scaled Video'); hold on;
    axis([1 scaledSize(2) 1 scaledSize(1)]); hold on;
    for i =1:stepSize:size(scaledAdjacencies,1)
        x=scaledPathlines(scaledAdjacencies(i,2),1,:); 
        y=scaledPathlines(scaledAdjacencies(i,2),2,:); 
        set(gca,'YDir','reverse');
        plot(x(:),y(:),'bx'); 
        hold on;
    end
    for i =1:stepSize:size(scaledAdjacencies,1)
        x=scaledPathlines(scaledAdjacencies(i,1),1,:); 
        y=scaledPathlines(scaledAdjacencies(i,1),2,:); 
        plot(x(:),y(:),'rx'); 
        hold on;
    end
    hold off;
    
    
end

%     PlotPathlines(originalPathlines,scaledPathlines,100,correspondances,...
%         originalAdjacencies,scaledAdjacencies , [vidHeight vidWidth] , newsize);


% function [correspondances] = GetCorrespondancesPerFrame(meshXWarped , meshYWarped ,...
%     pathlinesOriginal , pathlinesWarped , originalSize)
%     
%     actualIndicesOriginal = find(any(pathlinesOriginal,2));
%     actualIndicesWarped   = find(any(pathlinesWarped,2));
%     pathlinesOriginal     = pathlinesOriginal(any(pathlinesOriginal,2),:); %remove zero rows
%     pathlinesWarped       = pathlinesWarped(any(pathlinesWarped,2),:); %remove zero rows
%     
%     if( ~( isempty(pathlinesOriginal) || isempty(pathlinesWarped) ))
%         n = meshXWarped(end:end);
%         m = meshYWarped(end:end);
%         [ Y , X ] = size(meshXWarped);
% 
%         [ xRegularSmall, yRegularSmall] = meshgrid(double(1:n),double(1:m));
%         [ meshX , meshY ] = meshgrid( 1: floor(n/(X-1)) :n, ...
%                                       1: floor(m/(Y-1)) :m );
%         meshX(:,end) = n;
%         meshY(end,:) = m;
% 
%         newCoordsX = interp2( meshX, meshY, meshXWarped, xRegularSmall, yRegularSmall, 'spline');
%         newCoordsY = interp2( meshX, meshY, meshYWarped, xRegularSmall, yRegularSmall, 'spline');
% 
%         pathlinesOriginal(:,1) = floor( pathlinesOriginal(:,1) * ( n/originalSize(2)) );
%         pathlinesOriginal(:,2) = floor( pathlinesOriginal(:,2) * ( m/originalSize(1)) );
%         pathlinesOriginal(pathlinesOriginal<1) = 1;
%         
%         correspondances = [];
%         for i = 1: size(pathlinesOriginal,1)
%             warpedX = round(newCoordsX( pathlinesOriginal(i,2) , pathlinesOriginal(i,1) ));
%             warpedY = round(newCoordsY( pathlinesOriginal(i,2) , pathlinesOriginal(i,1) ));
%             warpedIndex = find( pathlinesWarped(:,1) < warpedX +1 & ...
%                                 pathlinesWarped(:,1) > warpedX -1 & ...
%                                 pathlinesWarped(:,2) < warpedY +1 & ...
%                                 pathlinesWarped(:,2) > warpedY -1 );
%             if ~isempty(warpedIndex)
%                 originalIndex = actualIndicesOriginal(i);
%                 warpedIndex = actualIndicesWarped( warpedIndex(1) );
%                 correspondances = [ correspondances; [originalIndex, warpedIndex ] ];
%             end
% 
%         end
%     end
% 
% end
