clc,clear

% Plot dispersion curve
moviePath = 'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00001.avi';
shotStart = 48;
shotEnd = 95;

video = VideoReader( moviePath );
frames = read(video);

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
%     plot(1:video.NumberOfFrames,dispersion);
plot(shotStart:shotEnd,dispersion,'LineWidth',1);
