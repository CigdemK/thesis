clear,clc;
F = dir('F:\Tez\Thesis\UserStudy\scaled2\*.mat');
for ii = 1:length(F)

%     video = VideoReader(['F:\Tez\Thesis\Hollywood2-actions\Hollywood2\UserStudy\' F(ii).name]);
%     frames = read(video);
%     frames = permute(frames.retargetedFrames,[4 2 3 1]);
    load(['F:\Tez\Thesis\UserStudy\scaled2\' F(ii).name]);
    [vidHeight, vidWidth, ~, nFrames] = size(retargetedFrames);
    for i = 1:nFrames
        a = im2uint8(retargetedFrames(:,:,:,i));
        newframes(:,:,:,i) = a;
%         frames(:,:,:,end+1) = cat(3,zeros(vidHeight,vidWidth),ones(vidHeight,vidWidth)*225,...
%             zeros(vidHeight,vidWidth));
    end
    
    writerObj = VideoWriter(['F:\Tez\Thesis\UserStudy\scaled2\' F(ii).name '.avi']);
    writerObj.FrameRate = 25;
    open(writerObj);
    writeVideo(writerObj,newframes);
    close(writerObj);

end

% ret = RetargetMethods;
% path = 'F:\Tez\Thesis\Hollywood2-actions\Hollywood2\SubsetAVIS\Frames\actioncliptest00186.avi_Frames.mat';
% load(path);
% retFrames = ret.Cigdem(path,frames,[16*25 9*25]);
