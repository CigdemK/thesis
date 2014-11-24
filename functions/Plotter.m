function F = Plotter
    F.PlotTrajectory2D = @PlotTrajectory2D;
    F.PlotTrajectory3D = @PlotTrajectory3D;
    F.PlotROC = @PlotROC;
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

function [score,R] = PlotROC(algorithm,plotColor)

    if nargin < 2
        plotColor = 'b';
    end
    if nargin == 0
        algorithm = 'sr';
    end

    Evaluate = EvalFunctions; 
    a = 1;
    for k = [397:597,607:824,862:924,947:1041,1072:1185,1213:1346,1502:1677,1678:1862,1871:1906,1937:2211,2246:2465,2469:2751]

        eyeFolder = strcat('../forLeMeur/massVideos/leMeurFixationMaps/frame_',num2str(k),'_saliencyMapGT.bmp');
        salFolder = strcat('../forLeMeur/massVideos/CalculatedSaliency_',algorithm,'/frame_',num2str(k),'.bmp');

        currentEye = imread(eyeFolder);
        currentSaliency = imread(salFolder);

        eyeFixation{a} = currentEye(:,:,1);
        saliency{a} = currentSaliency(:,:,1);
        a = a + 1;

    end

    [score,R] = Evaluate.CalculateAUCscoreVideo( saliency, eyeFixation, algorithm );

    plot(R(:,1),R(:,2),plotColor);

end
