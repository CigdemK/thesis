%% Init
clc,clear
moviePath = 'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00001.avi';
shotStart = 48;
shotEnd = 95;

video = VideoReader( moviePath );
frames = read(video);

%% Plot comparison line graph
RetargetingMethods = {'Rubinstein2008','Wang2011','Yan2013','Nei2013','Linear','Cigdem'};
ComparisonMethods = {'Compression Ratio','Waving','Blur','Sharpness','Jitter'};

nrRetMethods = size(RetargetingMethods,2);
nrCompMethods = size(ComparisonMethods,2);

Scores = rand([nrRetMethods nrCompMethods]).*100; % #Retmethods x #Compmethods

colorArr = {'b','g','k','r','m','c','y'}; 
figure;
for i = 1:nrRetMethods
    plot(Scores(i,:), 1:nrCompMethods,'Color',colorArr{mod(i,nrRetMethods)+1},'LineWidth',3);
    hold on;
end
set(gca, 'XTick', 0:20:100 , 'YTickLabel', ComparisonMethods);
title(['Comparison of Retaregeting Methods'],'FontWeight','Bold');
xlabel('Scores','FontWeight','Bold');
ylabel('Comparison Metrics','FontWeight','Bold');
legend(RetargetingMethods,'location','eastoutside');
%% Plot dispersion curve

dispersion = zeros(shotEnd-shotStart+1,1);
load('videoSaliency.mat');
tic;
saliencyMaps = mat2gray(videoSaliency);
for i = shotStart:shotEnd %1:video.NumberOfFrames

    currentSaliency = saliencyMaps(:,:,i-shotStart+1);
    [x,y] = find(currentSaliency > 0.5);
    nrOfPoints = size(x,1);
%     X(i-shotStart+1) = var(x);

    currentDispersion = 0;
    for k = 1: nrOfPoints
        for t = k+1:nrOfPoints
            currentDispersion = currentDispersion + (double((x(k)-x(t))^2 + (y(k)-y(t))^2));
        end
    end
    dispersion(i-shotStart+1) = currentDispersion / (nrOfPoints^2);
    toc;
end
    
figure;
plot(shotStart:shotEnd-1,dispersion(1:end-1),'LineWidth',2);
ylim([1 max(dispersion)]);
xlabel('Frames','FontWeight','Bold');
ylabel('Dispersion','FontWeight','Bold');
hold on;
p = polyfit(1:(shotEnd-shotStart),dispersion(1:end-1)',3);
y = polyval(p,1:shotEnd-shotStart);
plot(shotStart:shotEnd-1,y,'LineWidth',2,'Color','r');


