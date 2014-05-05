function F = UtilFunctions
    F.CheckIsInside = @CheckIsInside;
    F.PlotTrajectory2D = @PlotTrajectory2D;
    F.PlotTrajectory3D = @PlotTrajectory3D;
    F.CalculateMapping = @CalculateMapping; %wont work
    F.Homography = @Homography;
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

function PlotTrajectory2D( data , fitPower )
    x = double(data( : , 1 ));
    y = double(data( : , 2 ));
%     p = polyfit( x , y , fitPower );
%     plotx= linspace( min(x) , max(x));
%     hold on;
    c = jet(98);
    scatter(x,y,50,c);    
    hold on;
%     ploty = polyval( p , plotx );
%     plot(plotx, ploty, '-');
    axis([250 350 200 300]);
    hold on;
    
%     labels = {'10th frame','20th frame','30th frame','40th frame','10th frame','20th frame','30th frame','40th frame','10th frame','20th frame'};
    a = colorbar;
    set(gca, 'CLim', [1, 98]);
    set(a, 'YTick', 1:10:98); 
    set(a, 'YTickLabel', {'1st Frame','10th Frame','20th Frame','30th Frame','40th Frame','50th Frame','60th Frame','70th Frame','80th Frame','90th Frame'});
    xlabel('Frame Width');
    ylabel('Frame Height');
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

function [ H ] = Homography(img1,img2)
% Maps the points in the second image to the first image.
% Returns the homography matrix of img1-->img2

    I = single(rgb2gray(img1));
    J = single(rgb2gray(img2));
      
    % Apply SIFT and find matching points
    [features1,distances1] = vl_sift(I);
    [features2,distances2] = vl_sift(J);
    [matches ~] = vl_ubcmatch(distances1,distances2,5);

    % Find transformation matrix and apply transformation
    H = findHomography( ...
         [features1( 1 , matches( 1 , : ) ) ; features1( 2 , matches( 1 , : ) ) ],...
         [features2(1 , matches( 2 , : ));features2( 2 , matches( 2 , : ) ) ]);

end