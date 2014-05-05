function F = WANGFunctions
    F.ScaleStretch = @ScaleStretch;

end

function imgOut = ScaleStretch(imgIn, newsize)

%     % Compute significance map
%     Wb = saliency(imgIn);
%     [Gx,Gy] = imGradient(double(imgIn(:,:,1)));
%     Wa = sqrt(Gx.^2 + Gy.^2);
%     W = mat2gray(Wa.*Wb);

    % Construct the mesh grid
    [height width ~] = size(imgIn);
    [meshX,meshY] = meshgrid(1:width,1:height);
        
    Y = size(meshY,1);
    X = size(meshX,2);
    numberOfQuads = ( Y - 1 ) * ( X - 1 );  
    
    % Calculate grid weights
%     Wq = [];
%     for w = 1:X-1
%         for h = 1:Y-1
%             y = meshY(h,w):meshY(h+1,w);
%             x = meshX(h,w):meshX(h,w+1);
%             Wq = [ Wq ; mean( mean( W( y , x ) ) ) ];
%         end
%     end
%     Wq = Wq / (max(max(Wq)));
    
    % Start iteration for calculating V'
    meshXNew = ceil(meshX/2);
    meshYNew = ceil(meshY/2);
    
    imgOut = WarpImageWithMesh(imgIn , meshX , meshY , meshXNew , meshYNew);
    
    while ~CheckIfDone(meshXNew , meshYNew  , meshX , meshY )

        % Calculate S and l with the current estimates of V'
        S = zeros(numberOfQuads,1);
        l = zeros(numberOfQuads,4);
        for k = 1:numberOfQuads
            
            yInd = rem(k-1,Y-1)+1;
            xInd = floor((k-1)/(Y-1))+1;
            currentQuad = GetQuad(meshX, meshY, xInd, yInd);
            currentQuadNew = GetQuad(meshXNew, meshYNew, xInd, yInd);
                   
            top = 0;
            bottom = 0;
            for i = 1:4 %for each edge of a quad
                top = top +  ( currentQuad(mod(i-1,4)+1,:) - currentQuad(mod(i,4)+1,:) ) * ...
                    ( currentQuadNew(mod(i-1,4)+1,:) - currentQuadNew(mod(i,4)+1,:) )';
                bottom = bottom + sum( abs ( currentQuad(mod(i-1,4)+1,:) - currentQuad(mod(i,4)+1,:)) )^2;
            end
            S(k) = top / bottom;
            
            for i = 1:4 %for each edge of a quad
                l(k,i) =  sum ( abs( currentQuadNew(mod(i-1,4)+1,:) - currentQuadNew(mod(i,4)+1,:) ) ) / ...
                    sum( abs( currentQuad(mod(i-1,4)+1,:) - currentQuad(mod(i,4)+1,:) ) );
            end
        end

        F = @(X)FindTotalEnergy(X, meshX , meshY , S , l,Wq);
        x0 = cat( 3 , meshXNew , meshYNew );
        options = optimset('MaxIter', 1);
        mesh = fminunc( F , x0(:) , options ); 
        
        meshXNew = reshape(mesh( 1 : Y*X ) , [Y X] );
        meshYNew = reshape(mesh( Y*X+1 : end) ,[Y X] );
        
    end
   
    % Warp image with the new mesh
%     imgOut = WarpImageWithMesh(imgIn , meshX , meshY , meshXNew , meshYNew);
    
    
end

function D = FindTotalEnergy(meshNew , meshX , meshY ,S , l ,Wq) %meshNew is vectorized
    
    Y = size(meshY,1);
    X = size(meshX,2);
    meshXNew = reshape(meshNew( 1 : Y*X ) , [Y X] );
    meshYNew = reshape(meshNew( Y*X+1 : end) ,[Y X] );
    numberOfQuads = ( Y - 1 ) * ( X - 1 );

    Du = 0;
    Dt = 0;
    for k = 1:numberOfQuads
        
        yInd = rem(k-1,Y-1)+1;
        xInd = floor((k-1)/(Y-1))+1;
        currentQuad = GetQuad(meshX, meshY, xInd, yInd);
        currentQuadNew = GetQuad(meshXNew, meshYNew, xInd, yInd);

        DuCurrent = 0;
        DtCurrent = 0;
        for i = 1:4 %for each edge of a quad
            DuCurrent = DuCurrent + ( sum( abs( currentQuadNew(mod(i-1,4)+1,:) - currentQuadNew(mod(i,4)+1,:) )) - ...
                S(k) * sum( abs( currentQuad(mod(i-1,4)+1,:) - currentQuad(mod(i,4)+1,:) ) ) )^2;
            DtCurrent = DtCurrent + ( sum( abs( currentQuadNew(mod(i-1,4)+1,:) - currentQuadNew(mod(i,4)+1,:) )) - ...
                l(k,i) * sum( abs( currentQuad(mod(i-1,4)+1,:) - currentQuad(mod(i,4)+1,:) ) ) )^2;
        end

        Du = Du + DuCurrent * Wq(k);
        Dt = Dt + DtCurrent * Wq(k);
        D = Du + Dt;
    end
    
end

function result = CheckIfDone( meshXNew , meshYNew  , meshX , meshY )

    X = abs(meshXNew - meshX);
    Y = abs(meshYNew - meshY);
    result = true;
    
    if size(X(X>0.5)) + size(Y(Y>0.5)) > 0
        result = false;
    end
    
end

function quad = GetQuad( meshX, meshY, xInd, yInd )
    quad = [[ meshY( yInd   , xInd ) meshX( yInd , xInd   )] ; ...
            [ meshY( yInd+1 , xInd ) meshX( yInd , xInd   )] ; ...
            [ meshY( yInd   , xInd ) meshX( yInd , xInd+1 )] ; ...
            [ meshY( yInd+1 , xInd ) meshX( yInd , xInd+1 )]];
end

function warpedImage = WarpImageWithMesh(imgIn , meshX , meshY , meshXNew , meshYNew)

    Image = double(imgIn);

    [m, n , ~] = size(Image);

    xi = 1:n;
    yi = 1:m;

    % calls the helper function getPoints to vectors of the 
    [x y zx, zy] = getPoints(meshX , meshY , meshXNew , meshYNew);

    % interpolate all the points in the image based on the displacements calculated above
    displaceX = round(griddata(x, y, zx, xi', yi, 'cubic'));   
    displaceY = round(griddata(x, y, zy, xi', yi, 'cubic'));

    % removes and NaN from the griddata interpolation
    displaceX(isnan(displaceX)) = 0; 
    displaceY(isnan(displaceY)) = 0; 

    % change back into matricies
    displaceX = reshape(displaceX, m,n);
    displaceY = reshape(displaceY, m,n);

    % calc where each pixel should be remapped
    [coordsX coordsY] = meshgrid(1:n, 1:m);

    newCoordsX = coordsX + displaceX;
    newCoordsY = coordsY + displaceY;

    % adjust any pixels that might be out of range of the image
    newCoordsX(newCoordsX < 1) = 1; 
    newCoordsX(newCoordsX > n) = n; 
    newCoordsY(newCoordsY < 1) = 1; 
    newCoordsY(newCoordsY > m) = m; 

    linearX = reshape(newCoordsX, 1, m*n); 
    linearY = reshape(newCoordsY, 1, m*n); 

    indecies = sub2ind([m,n], linearY, linearX);  
    Image = reshape(Image,1, n*m,3);
    warpedImage = Image(1,indecies,:);
    warpedImage = reshape(uint8(warpedImage), m, n,3);   % reshape back into an image

end

function [meshX, meshY, zx, zy]   =  getPoints(meshX , meshY , meshXNew , meshYNew)

    Y = size(meshY,1);
    X = size(meshX,2);
    newSize = [ meshX(Y,X) , meshY(Y,X) ];
    
    meshX = meshX(:);
    meshY = meshY(:);
    meshXNew = meshXNew(:);
    meshYNew = meshYNew(:);
    
    numVertices = Y * X;  
    
    zx = zeros(numVertices,1);
    zy = zeros(numVertices,1);
    for i = 1:numVertices
        zx(i) = (meshX(i) - meshXNew(i))*(meshX(i)/meshXNew(i));   % x offset component of the displacement vector
        zy(i) = (meshY(i) - meshYNew(i))*(meshY(i)/meshYNew(i));   % x offset component of the displacement vector
    end

    % add corners before for better interpolation
%     meshX(i+1) = 1;             meshY(i+1) = 1;             zx(i+1) = 0;    zy(i+1) = 0;
%     meshX(i+2) = 1;             meshY(i+2) = newSize(2);    zx(i+2) = 0;    zy(i+2) = 0;        
%     meshX(i+3) = newSize(1);    meshY(i+3) = 1;             zx(i+3) = 0;    zy(i+3) = 0;
%     meshX(i+4) = newSize(1);    meshY(i+4) = newSize(2);    zx(i+4) = 0;    zy(i+4) = 0;
    
end



