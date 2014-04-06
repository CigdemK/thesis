function F = UtilFunctions
    F.CheckIsInside = @CheckIsInside;
    F.PlotTrajectory2D = @PlotTrajectory2D;
    F.PlotTrajectory3D = @PlotTrajectory3D;
    F.CalculateMapping = @CalculateMapping; %wont work
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
