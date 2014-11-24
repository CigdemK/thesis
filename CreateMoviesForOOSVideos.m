clc, clear;

moviePath = 'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptrain00438.avi';
regexres = regexp(moviePath,'.avi','split');
myPath = regexres{1};

% Read movie
video = VideoReader( moviePath );
originalFrames = read(video);

% Frame definitions
OOSFrames = zeros(video.Height,video.Width,3,video.NumberOfFrames);
normalFrames = zeros(video.Height,video.Width,3,video.NumberOfFrames);
frames = zeros(video.Height*3+100,video.Width,3,video.NumberOfFrames);

for i = 1:video.NumberOfFrames

    % Read normal and oos heatmaps
    normalFrames(:,:,:,i) = im2double(imread([myPath '\LeMeurNormalResults\frame_' num2str(i) '_heatMap.bmp']));
    OOSFrames(:,:,:,i) = im2double(imread([myPath '\LeMeurOOSResults\frame_' num2str(i) '_heatMap.bmp']));
   
    % Read normal fixations from .stat files & parse coordinates
    normalFixations = fileread([myPath '\LeMeurNormalStats\frame_' num2str(i) '.stat']);
    normalFixationsX = []; normalFixationsY = [];
    reg = regexp(normalFixations,'\n','split');
    for k = 1:size(reg,2)
        regex = regexp(reg{k},' ','split');
        normalFixationsX = [normalFixationsX , regex(1:3:end-3)];
        normalFixationsY = [normalFixationsY , regex(2:3:end-3)];
    end
    normalFixationsX = normalFixationsX(:); normalFixationsY = normalFixationsY(:);
    
    % Read oos fixations from .stat files & parse coordinates
    OOSFixations = fileread([myPath '\LeMeurOOSStats\frame_' num2str(i) '.stat']);
    OOSFixationsX = [] ; OOSFixationsY = [];
    reg = regexp(OOSFixations,'\n','split');
    for k = 1:size(reg,2)
        regex = regexp(reg{k},' ','split');
        OOSFixationsX = [OOSFixationsX , regex(1:3:end-3)];
        OOSFixationsY = [OOSFixationsY , regex(2:3:end-3)];
    end
    OOSFixationsX = OOSFixationsX(:);  OOSFixationsY = OOSFixationsY(:);
    
    % Mark normal fixations with white dot
    for k = 1:size(normalFixationsX)
        originalFrames(str2num(normalFixationsY{k}):str2num(normalFixationsY{k})+4 , ...
            str2num(normalFixationsX{k}):str2num(normalFixationsX{k})+4  ,:, i) = 254*ones() ;   
    end
    
    % Mark oos fixation with red dot
    redDot = cat(3, 254*ones(5,5), zeros(5,5), zeros(5,5));
    for k = 1:size(OOSFixationsX)
        originalFrames(str2num(OOSFixationsY{k}):str2num(OOSFixationsY{k})+4 , ...
            str2num(OOSFixationsX{k}):str2num(OOSFixationsX{k})+4  ,:, i) = redDot;   
    end
    originalFrames = originalFrames(1:video.Height,1:video.Width,:,:);
    
    % Vertically concatenate all 3 frames
    frames(:,:,1,i) = [im2double(originalFrames(:,:,1,i)); zeros(50,video.Width); normalFrames(:,:,1,i); zeros(50,video.Width); OOSFrames(:,:,1,i)];
    frames(:,:,2,i) = [im2double(originalFrames(:,:,2,i)); zeros(50,video.Width); normalFrames(:,:,2,i); zeros(50,video.Width); OOSFrames(:,:,2,i)];
    frames(:,:,3,i) = [im2double(originalFrames(:,:,3,i)); zeros(50,video.Width); normalFrames(:,:,3,i); zeros(50,video.Width); OOSFrames(:,:,3,i)];
end


writerObj = VideoWriter([myPath '\Results.avi']);
writerObj.FrameRate = video.FrameRate;
open(writerObj);
writeVideo(writerObj,frames);
close(writerObj);