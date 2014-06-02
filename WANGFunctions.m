function F = WANGFunctions
    F.ScaleStretch = @ScaleStretch;
    F.PerFrameOptimization = @PerFrameOptimization;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Public Functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scaledPathlines = PerFrameOptimization(FolderName, FileName, newsize)

%     TimerScript;
%     
% %     Read input movie
%     Reader = ReadFunctions;
%     Cropper = CropFunctions;
%     [video,nFrames,vidHeight,vidWidth,vidFPS,frames] = Reader.ReadData(FolderName,FileName);
%     mov = Reader.NewMovie(nFrames, vidHeight,vidWidth);
%     mov = Reader.ReadMovie(mov,video);
%     
% %     Remove black bars and adjust the size parameters
%     movNoBars = Cropper.RemoveBlackBars(mov);
%     [vidHeight,vidWidth,~] = size(movNoBars(1).cdata);
%     
% % Step 1: Apply Scale&Stretch to all frames
%     tic;
%     originalFrames = zeros(vidHeight,vidWidth,3,nFrames);
%     scaledFrames = zeros(newsize(1),newsize(2),3,nFrames);
%     gridX = zeros(ceil(vidHeight/50),ceil(vidWidth/50),nFrames);
%     gridY = zeros(ceil(vidHeight/50),ceil(vidWidth/50),nFrames);
%     for k = 1:nFrames
%         originalFrames(:,:,:,k) = im2double(movNoBars(k).cdata);
%         [tmp, gridX(:,:,k), gridY(:,:,k)]= ScaleStretch(movNoBars(k).cdata,newsize);
%         scaledFrames(:,:,:,k) = im2double(tmp);
%     end
%     toc;
%     
% % Step 2: Calculate pathlines, adjacency and correspondance
% 
%     originalPathlines   = GetPathlines( originalFrames, [vidHeight/20 vidWidth/20] ); toc;
%     scaledPathlines     = GetPathlines( scaledFrames, [vidHeight/20 vidWidth/20] ); toc;
% 
%     originalAdjacencies = CalculateAdjacencies(originalPathlines); toc;
%     scaledAdjacencies   = CalculateAdjacencies(scaledPathlines);  toc;
% 
%     correspondances = GetCorrespondances(gridX , gridY , originalPathlines , scaledPathlines, [vidHeight, vidWidth] );
  
%     save('tmp.mat'); 
    load('tmp.mat');
    
% % Step 3: Optimize Scaling with the Pathlines   
    
    
    
%     PlotPathlines(originalPathlines,scaledPathlines,100,correspondances,...
%         originalAdjacencies,scaledAdjacencies , [vidHeight vidWidth] , newsize);
end

function [imgOut, meshXNew , meshYNew]= ScaleStretch(imgIn, newsize)

    % Compute significance map
    Wb = gbvs(imgIn);
    [Gx,Gy] = imGradient(double(imgIn(:,:,1)));
    Wa = sqrt(Gx.^2 + Gy.^2);
    W = mat2gray(Wa.*Wb.master_map_resized);

    % Construct the mesh grid
    [height width ~] = size(imgIn);
    [meshX,meshY] = meshgrid(1:50:width,1:50:height);
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
    % matrices are procomputed to be used in all iterations.
    mask = VertexMask(X,Y);
    originalVetricesX = mask*meshX(:);
    originalVerticesY = mask*meshY(:);

    while ~CheckIfDone(meshXOld , meshYOld  , meshXNew , meshYNew )

        meshXOld = reshape(meshXNew, [Y X] );
        meshYOld = reshape(meshYNew, [Y X] );
        
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

        FX = @(Vnew)sum(((mask*Vnew).^2) .* (Wq + 1) + ...
                         ((mask*meshX(:)).^2) .* (Wq .* (S.^2) + (l.^2)) + ...
                         (mask*Vnew).*(mask*meshX(:)) .* ((Wq .* S + l) * (-2)));
        FY = @(Vnew)sum(((mask*Vnew).^2) .* (Wq + 1) + ...
            ((mask*meshY(:)).^2) .* (Wq .* (S.^2) + (l.^2)) + ...
            (mask*Vnew).*(mask*meshY(:)) .* ((Wq .* S + l) * (-2)));
        
        % Set limits of fmincon so that the output image will be
        % rectangular
        options = optimset('Algorithm','interior-point','Display','off');
        AeqX = zeros(Y*X);
        beqX = zeros(Y*X,1);
        AeqY = zeros(Y*X);
        beqY = zeros(Y*X,1);
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
        meshXNew = fmincon( FX ,meshXOld(:) ,[] ,[] ,AeqX ,beqX ,[] ,[] ,[], options );
        meshYNew = fmincon( FY ,meshYOld(:) ,[] ,[] ,AeqY ,beqY ,[] ,[] ,[], options );

    end

    % Warp image with the new mesh
    meshXNew = reshape(meshXNew,Y,X);
    meshYNew = reshape(meshYNew,Y,X);
    imgOut = WarpImageWithMesh(imgIn , meshXNew , meshYNew);    
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Private Functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Functions for PerFrameOptimization

function [pathlines] = GetPathlines(frames , meshSize)
    
    [vidHeight , vidWidth, ~, nFrames] = size(frames);
    
    % Calculate optical flow
    [Vx, Vy] = CalculateOpticalFlow(frames);
%     load('tmp2_gooddouble.mat');
%     save('tmp2_gooddouble.mat','Vx','Vy');
    
    % Initialize mesh grids
    [ initMeshX , initMeshY ] = meshgrid( 1:floor(vidWidth/meshSize(2)):vidWidth, ...
                                          1:floor(vidHeight/meshSize(1)):vidHeight );
                                      
    X = size(initMeshX,2);
    Y = size(initMeshX,1);
    
    meshX = zeros(Y,X,nFrames);
    meshY = zeros(Y,X,nFrames);
    
    meshX(:,:,1) = initMeshX;
    meshY(:,:,1) = initMeshY;

    % Calculate pathlines
    maxNrOfPathlines = Y*X*nFrames;
    pathlines = zeros(maxNrOfPathlines, 2, nFrames);
    [ xReg , yReg ] = meshgrid(1:vidWidth,1:vidHeight);
    j = 1;
    for k = 2:nFrames
               
        % Add the current optical flow to the prev. vertex pts
        xPts = interp2(xReg , yReg , Vx(:,:,k-1) , meshX(:,:,k-1) , meshY(:,:,k-1) );
        yPts = interp2(xReg , yReg , Vy(:,:,k-1) , meshX(:,:,k-1) , meshY(:,:,k-1) );
        meshX(:,:,k) = meshX(:,:,k-1) + xPts;
        meshY(:,:,k) = meshY(:,:,k-1) + yPts;
        
        % Find terminated paths that exceeds frame size
        [indXr indXc] = find( meshX(:,:,k)<1 | meshX(:,:,k)>vidWidth );
        [indYr indYc] = find( meshY(:,:,k)<1 | meshY(:,:,k)>vidHeight );
        
        allRows     = [indXr;indYr];
        allColumns  = [indXc;indYc];

        equals=[];
        for i = 1:size(allRows)
            equalRows = find( allRows(i) == allRows(i+1:end) ) + i;
            equalCols = find( allColumns(i) == allColumns(i+1:end) ) + i;
            equals = [equals;intersect(equalRows , equalCols)];

        end
        allRows(equals) = [];
        allColumns(equals) = [];

        % Add the terminated paths to path list and feed new paths from
        % initMesh
        if( ~(isempty(allRows)) )
               
            % add them to pahtlines
            for i = 1:size(allRows)
                
                meshX(allRows(i),allColumns(i),k) = 0;
                meshY(allRows(i),allColumns(i),k) = 0;
                
                if size(find(meshX(allRows(i),allColumns(i),:) > 0),1) > 2 % dont add the pathline if there are less than 3 point in it
                    
                    pathlines(j,1,:) = meshX(allRows(i),allColumns(i),:);
                    pathlines(j,2,:) = meshY(allRows(i),allColumns(i),:);
                    j = j + 1;
                    
                end
                
                meshX(allRows(i),allColumns(i),:) = 0;
                meshY(allRows(i),allColumns(i),:) = 0;
                meshX(allRows(i),allColumns(i),k) = initMeshX(allRows(i),allColumns(i));
                meshY(allRows(i),allColumns(i),k) = initMeshY(allRows(i),allColumns(i));
                
            end
            
        end
        
    end
    
    for i = 1:Y
        for k = 1:X
            if size(find(meshX(i,k,:) > 0),1) > 2 
                pathlines(j,1,:) = meshX(i,k,:);
                pathlines(j,2,:) = meshY(i,k,:);
                j = j + 1;
            end
        end
    end
    
    % remove all zero lines
    pathlines(all(all(pathlines == 0,3),2),:,:) = [];
    
end

function [Vx, Vy] = CalculateOpticalFlow(frames)
    
    [vidHeight , vidWidth, ~, nFrames] = size(frames);
    
    Vx = zeros(vidHeight,vidWidth,nFrames);
    Vy = zeros(vidHeight,vidWidth,nFrames);
    
    for k = 1:nFrames-1
        [VxTmp, VyTmp, ~] = Coarse2FineTwoFrames(frames(:,:,:,k),frames(:,:,:,k+1));
        Vx(:,:,k) = VxTmp;
        Vy(:,:,k) = VyTmp;
    end
    
    Vx(:,:,end) = Vx(:,:,end-1);
    Vy(:,:,end) = Vy(:,:,end-1);
    
    Vx = Vx.*50;
    Vy = Vy.*50;
    
end

function [adjacencies] = CalculateAdjacencies(pathlines)

    adjacencies = [];
    [ ~, ~, nFrames ] = size(pathlines);
    
    for i = 1:nFrames
        
        timeSlice = pathlines(:,:,i);
        actualIndices = find(any(timeSlice,2));
        timeSlice = timeSlice(any(timeSlice,2),:); %remove zero rows

        if(i~=1)
            prevTimeSlice = pathlines(:,:,i-1);
            prevIndices = find(any(prevTimeSlice == 0 ,2)); 
            [ actualIndices , ~ , tmp] = intersect( prevIndices , actualIndices);
            timeSlice = timeSlice(tmp , :);
        end

        for k = 1:size(timeSlice,1)
            
            X = k+find( abs( timeSlice( k+1:end, 1 ) - timeSlice( k, 1) ) < 21 );
            Y = k+find( abs( timeSlice( k+1:end, 2 ) - timeSlice( k, 2) ) < 21 );
            currentAdj = intersect(X,Y);
            
            actualIndexSecond = actualIndices(currentAdj);
            actualIndexFirst = repmat( actualIndices(k) , size(currentAdj) );
            adjacencies = [adjacencies; [ actualIndexFirst actualIndexSecond repmat(i,size(currentAdj))] ] ;
        
        end
    end
    
end

function [correspondances] = GetCorrespondances(meshXWarped , meshYWarped , ...
    originalPathlines , scaledPathlines , originalSize)

    correspondances = [];
    [~,~,nFrames] = size(meshXWarped);
    
    for k = 1:nFrames
        
        timeSliceOriginal = originalPathlines(:,:,k);
        timeSliceScaled = scaledPathlines(:,:,k);
        indOriginal = 1:size(timeSliceOriginal,1);
        indScaled = 1:size(timeSliceScaled,1);
        
        currentCorr = GetCorrespondancesPerFrame(meshXWarped(:,:,k) , meshYWarped(:,:,k),...
            timeSliceOriginal , timeSliceScaled , originalSize);
        correspondances =  [correspondances; currentCorr ];
        
    end
    
    correspondances = unique(correspondances, 'rows');
    
end

function [correspondances] = GetCorrespondancesPerFrame(meshXWarped , meshYWarped ,...
    pathlinesOriginal , pathlinesWarped , originalSize)
    
    actualIndicesOriginal = find(any(pathlinesOriginal,2));
    actualIndicesWarped   = find(any(pathlinesWarped,2));
    pathlinesOriginal     = pathlinesOriginal(any(pathlinesOriginal,2),:); %remove zero rows
    pathlinesWarped       = pathlinesWarped(any(pathlinesWarped,2),:); %remove zero rows
    
    if( ~( isempty(pathlinesOriginal) || isempty(pathlinesWarped) ))
        n = int32(meshXWarped(end:end));
        m = int32(meshYWarped(end:end));
        [ Y , X ] = size(meshXWarped);

        [ xRegularSmall, yRegularSmall] = meshgrid(double(1:n),double(1:m));
        [ meshX , meshY ] = meshgrid( 1 : ceil( double(n/int32(X)) ) : double(n) , ...
                                      1 : ceil( double(m/int32(Y)) ) : double(m) );
        meshX(:,end) = double(n);
        meshY(end,:) = double(m);

        newCoordsX = interp2( meshX, meshY, meshXWarped, xRegularSmall, yRegularSmall, 'cubic');
        newCoordsY = interp2( meshX, meshY, meshYWarped, xRegularSmall, yRegularSmall, 'cubic');

        pathlinesOriginal(:,1) = ceil( pathlinesOriginal(:,1) * ( double(n) / originalSize(2) ) );
        pathlinesOriginal(:,2) = ceil( pathlinesOriginal(:,2) * ( double(m) / originalSize(1) ) );

        correspondances = [];
        for i = 1: size(pathlinesOriginal,1)
            warpedX = round(newCoordsX( pathlinesOriginal(i,2) , pathlinesOriginal(i,1) ));
            warpedY = round(newCoordsY( pathlinesOriginal(i,2) , pathlinesOriginal(i,1) ));
            warpedIndex = find( pathlinesWarped(:,1) < warpedX +1 & ...
                                pathlinesWarped(:,1) > warpedX -1 & ...
                                pathlinesWarped(:,2) < warpedY +1 & ...
                                pathlinesWarped(:,2) > warpedY -1 );
            if ~isempty(warpedIndex)
                originalIndex = actualIndicesOriginal(i);
                warpedIndex = actualIndicesWarped( warpedIndex(1) );
                correspondances = [ correspondances; [originalIndex, warpedIndex ] ];
            end

        end
    end

end

% % Functions for ScaleStretch

function binaryMatrix =  VertexMask(X,Y)

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

function result = CheckIfDone( meshXOld , meshYOld  , meshX , meshY )

    X = abs(meshXOld(:) - meshX(:));
    Y = abs(meshYOld(:) - meshY(:));
    result = true;
    
    if size(X(X>0.5)) + size(Y(Y>0.5)) > 0
        result = false;
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
    n = int32(meshX(end:end));
    m = int32(meshY(end:end));
    
    imgIn = imresize(imgIn,[m,n]);

    [xRegular,yRegular] = meshgrid(double(1:n),double(1:m));
    [meshXRegular,meshYRegular] = meshgrid( 1 : ceil( double(n/int32(X)) ) : double(n) , ...
                                            1 : ceil( double(m/int32(Y)) ) : double(m) );
    meshXRegular(:,end) = double(n);
    meshYRegular(end,:) = double(m);
    
    newCoordsX = interp2( meshXRegular, meshYRegular, meshX, xRegular, yRegular, 'cubic');
    newCoordsY = interp2( meshXRegular, meshYRegular, meshY, xRegular, yRegular, 'cubic');
    
    warped1 = interp2( xRegular, yRegular, imgIn(:,:,1), newCoordsX, newCoordsY, 'cubic');
    warped2 = interp2( xRegular, yRegular, imgIn(:,:,2), newCoordsX, newCoordsY, 'cubic');
    warped3 = interp2( xRegular, yRegular, imgIn(:,:,3), newCoordsX, newCoordsY, 'cubic');
    
    warpedImage = uint8(cat(3,warped1,warped2,warped3));

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



