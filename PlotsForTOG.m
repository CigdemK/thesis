% % %% Init
% % % clc,clear
% % moviePath = 'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00001.avi';
% % shotStart = 48;
% % shotEnd = 95;
% % 
% % video = VideoReader( moviePath );
% % frames = read(video);
% 
% %% Plot comparison line graph
% % RetargetingMethods = {'Rubinstein2008','Yan2013','Linear','Cigdem','Original'};
% % ComparisonMethods = {'Shakiness','Focus','Blur','Sharpness','Brightness','Compress','PSNR','MSE', ...
% %                 'VIFP','VSNR','MSSIM','SSIM','UQI','SNR','WSNR','MDE'};
% % 
% % RetargetingMethods = {'Original'};
% % ComparisonMethods = {'UQI','Blur','Focus','Sharpness'};
% % 
% % nrRetMethods = size(RetargetingMethods,2);
% % nrCompMethods = size(ComparisonMethods,2);
% % OriginalScores = [0.4806 8.3533e+17 0.4288 773.9010 0.1291 60.4559 27.9595 155.1156 0.3998 19.9915 0.8876 0.8411 0.5667 13.4895 16.6465 0.5792];
% % Scores = [1.5639 8.9021e+17 0.4251 699.3531 0.1368 48.4364 23.2700 391.0211 0.2562 12.5505 0.7506 0.7206 0.3718 9.1319 10.9799 0.2220;
% %           0.6244 9.0215e+17 0.4277 618.8033 0.1383 49.7018 26.9505 185.0955 0.3527 18.9200 0.8711 0.8106 0.5336 12.8754 15.8865 0.1889;
% %           0.2726 8.4255e+17 0.4192 708.0024 0.1291 51.0144 27.9772 154.3458 0.3886 20.2987 0.8885 0.8345 0.5549 13.5044 16.6365 0.6370;
% %           1.1120 5.6954e+17 0.4284 781.9192 0.1133 59.3698 29.8146 163.3683 0.3496 19.2704 0.8855 0.8791 0.3876 13.4617 15.3687 0.7942;
% %           0.4806 8.3533e+17 0.4288 773.9010 0.1291 60.4559 27.9595 155.1156 0.3998 19.9915 0.8876 0.8411 0.5667 13.4895 16.6465 0.5792]; % #Retmethods x #Compmethods
% % 
% %       Scores = [8.9021e+17 0.4251 699.3531 0.1368 48.4364 23.2700 391.0211 0.2562 12.5505 0.7506 0.7206 0.3718 9.1319 10.9799;
% %           9.0215e+17 0.4277 618.8033 0.1383 49.7018 26.9505 185.0955 0.3527 18.9200 0.8711 0.8106 0.5336 12.8754 15.8865;
% %           8.4255e+17 0.4192 708.0024 0.1291 51.0144 27.9772 154.3458 0.3886 20.2987 0.8885 0.8345 0.5549 13.5044 16.6365;
% %           1.3768e+18 1.1011 315.1053 0.2581 44.8063 45.9941 816.1684 0.4182 20.1859 1.3963 1.4567 0.5226 17.3710 19.3270]; % #Retmethods x #Compmethods
% % Scores(1,:) = Scores(1,:) ./ OriginalScores;
% % Scores(2,:) = Scores(2,:) ./ OriginalScores;
% % Scores(3,:) = Scores(3,:) ./ OriginalScores;
% % Scores(4,:) = Scores(4,:) ./ OriginalScores;
% % Scores(5,:) = Scores(5,:) ./ OriginalScores;      
% % 
% % 
% colorArr = {'b','g','k','r','m','c','y'}; 
% figure;
load('data\D.mat');
ComparisonMethods = D(1).methods;
nrCompMethods = size(ComparisonMethods,2);
maxarray = zeros(1,nrCompMethods);
minarray = ones(1,nrCompMethods)*10^30;
maxOutliers = zeros(1,nrCompMethods);
minOutliers = zeros(1,nrCompMethods);
for j = 1:nrCompMethods
    for i = 1:230
        if isinf(D(i).results(j)) 
            D(i).results(j) = 0;
        end
        if(D(i).results(j)>maxarray(j))
            maxarray(j)=D(i).results(j);
            maxOutliers(j)=i;
        end
        if(D(i).results(j)<minarray(j))
            minarray(j)=D(i).results(j);
            minOutliers(j)=i;
        end
    end
end
% % for i = 1:nrRetMethods
%     for j = 1:size(D,2)
%         currentRes = D(j).results;
%         currentRes(isnan(currentRes)) = 0;
%         currentRes = currentRes./maxarray;
%         plot(currentRes, 1:nrCompMethods,'Color',colorArr{mod(j,7)+1},'LineWidth',1);
% %         plot(Scores(i,:), 1:nrCompMethods,'Color',colorArr{mod(i,nrRetMethods)+1},'LineWidth',2);
%         hold on;
%     end
% % end
% 
% set(gca,'YTick', 1:nrCompMethods ,'XTick', 1:100 , 'YTickLabel', ComparisonMethods);
% % set(gca, 'XTick', 0:20:100 , 'YTickLabel', ComparisonMethods);
% title('Comparison of Retaregeting Methods','FontWeight','Bold');
% xlabel('Scores','FontWeight','Bold');
% ylabel('Comparison Metrics','FontWeight','Bold');
% % legend(RetargetingMethods,'location','eastoutside');
% axis tight;
% 
% 
% %% Plot dispersion curve
% % 
% % dispersion = zeros(shotEnd-shotStart+1,1);
% % load('videoSaliency.mat');
% % tic;
% % saliencyMaps = mat2gray(videoSaliency);
% % for i = shotStart:shotEnd %1:video.NumberOfFrames
% % 
% %     currentSaliency = saliencyMaps(:,:,i-shotStart+1);
% %     [x,y] = find(currentSaliency > 0.5);
% %     nrOfPoints = size(x,1);
% % %     X(i-shotStart+1) = var(x);
% % 
% %     currentDispersion = 0;
% %     for k = 1: nrOfPoints
% %         for t = k+1:nrOfPoints
% %             currentDispersion = currentDispersion + (double((x(k)-x(t))^2 + (y(k)-y(t))^2));
% %         end
% %     end
% %     dispersion(i-shotStart+1) = currentDispersion / (nrOfPoints^2);
% %     toc;
% % end
% % shotStart = 48;
% % shotEnd = 95;
% % load('dispersion.mat')
% % h1 = figure;
% % plot(shotStart:shotEnd,dispersion(1:end-1),'LineWidth',2);
% % ylim([1 max(dispersion)]);
% % xlabel('Frames','FontWeight','Bold');
% % ylabel('Dispersion','FontWeight','Bold');
% % hold on;
% % p = polyfit(1:(shotEnd-shotStart)+1,dispersion(1:end-1)',3);
% % y = polyval(p,1:shotEnd-shotStart+1);
% % plot(shotStart:shotEnd,y,'LineWidth',2,'Color','r');
% % print(h1,'-djpeg','-r150','dispersionGraph.jpg');
% 
% % %% Plot video cube with trajectories
% % clc,clear;
% % folderName = 'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00001.avi';
% % 
% % video=VideoReader(folderName);
% % frames = read(video);
% % [height,width,~,nFrames] = size(frames);
% % load ([folderName '_impTrajectories.mat'])
% % shotStart = 48;
% % shotEnd = 95;
% % 
% % % Plot video cube
% % h1 = figure;
% % for i = shotStart:shotEnd
% % %     img = imread(['F:\Thesis\bin\data\actioncliptest00001\LeMeurNormalResults\frame_' num2str(i) '_heatMap.bmp']);
% %     img = frames(:,:,:,i);
% %     cdata = flipdim( img, 1 );
% %     surface([0 width; 0 width], [(i-1)*10 (i-1)*10; (i-1)*10 (i-1)*10], [0 0; height height], ...
% %         'FaceColor', 'texturemap', 'CData', cdata , 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% %     set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
% %     xlabel('Video Width','FontWeight','Bold');
% %     ylabel('Time / Frames','FontWeight','Bold');
% %     zlabel('Video Height','FontWeight','Bold');
% %     axis equal;
% %     view([-34 68])
% % end
% % 
% % % Plot all trajectories
% % for k = 1:100:size(impTrajectories,2)
% % 
% %     currentTrajectory = impTrajectories(k);
% %     tr = currentTrajectory.trajectory;
% %    
% %     trajectoryLength = size(tr,2);
% % 
% %     trajectoryEnd = currentTrajectory.frameNum;
% %     trajectoryStart = trajectoryEnd - trajectoryLength + 1;
% %     if (trajectoryEnd > shotEnd) || (trajectoryStart < shotStart); continue; end
% %     
% %     hold on;
% %     plot3(tr(2,:),(trajectoryStart-1)*10:10:(trajectoryEnd-1)*10,tr(1,:),'Color','r','linewidth',3);
% % %     plot(tr(2,:),(trajectoryStart-1)*10:10:(trajectoryEnd-1)*10,'x')
% % end
% % 
% % print(h1,'-djpeg','-r150','test.jpg');










