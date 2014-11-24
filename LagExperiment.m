clc,clear

folderName = 'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\';
AITE = ActionsInTheEye;
shotStart = 48;
shotEnd = 95;
    
for ultimatenumber = 1:1
    
    fileName = ['actioncliptest0000' num2str(ultimatenumber)  '.avi'];
    moviePath = [folderName fileName];
    
    video = VideoReader( moviePath );
    frames = read(video);
    frames = frames(:,:,:,shotStart:shotEnd);
    
    resultMap = AITE.ReadEyeTrackingData(moviePath);
    resultMap.vidDuration = video.Duration;
    resultMap.nFrames = video.NumberOfFrames;
    saliencyPoints = AITE.CalculateMapping(resultMap);

    dispersion = zeros(video.NumberOfFrames,1);
    load('saliencyMaps.mat');
    tic;
    saliencyMaps = mat2gray(saliencyMaps);
    for i = shotStart:shotEnd %1:video.NumberOfFrames

        indices = find(saliencyPoints( : , 2 ) == i);
        x = double(saliencyPoints(indices,4));
        y = double(saliencyPoints(indices,3));

%         currentSaliency = saliencyMaps(:,:,i-shotStart+1);
%         X = find(currentSaliency > 0);
        nrOfPoints = size(x,1);
%         X(i-shotStart+1) = var(x);
        currentVariance = 0;
        for k = 1: nrOfPoints
            for t = k+1:nrOfPoints
                currentVariance = currentVariance + (double((x(k)-x(t))^2 + (y(k)-y(t))^2));
            end
        end
        dispersion(i) = currentVariance / (nrOfPoints^2);
        toc;
    end

    figure;
%     plot(1:video.NumberOfFrames,dispersion);
    plot(shotStart:shotEnd,dispersion(shotStart:shotEnd));
%     title(fileName);
%     ylabel('Dispersion');
%     xlabel('Frame Number');
    
end
    
