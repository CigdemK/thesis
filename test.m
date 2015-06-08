clear,clc;
F = dir('E:\Tez\Thesis\userStudy\originalFrames\nguyen\*.mat');
for ii = 1:length(F)

%     video = VideoReader(['E:\Tez\Thesis\userStudy\originalAvi\' F(ii).name]);
%     frames = read(video);
%     frames = permute(frames.retargetedFrames,[4 2 3 1]);
    frames2 = [];
    load(['E:\Tez\Thesis\userStudy\originalFrames\nguyen\' F(ii).name]);
    [vidHeight, vidWidth, ~, nFrames] = size(frames2);

    for i = 1:25
%         a = im2uint8(retargetedFrames(:,:,:,i));
%         newframes(:,:,:,i) = a;
        frames2(:,:,:,end+1) = cat(3,zeros(vidHeight,vidWidth),ones(vidHeight,vidWidth)*225,...
            zeros(vidHeight,vidWidth));
    end
    
    writerObj = VideoWriter(['E:\Tez\Thesis\userStudy\originalFrames\nguyen\' F(ii).name '.avi']);
    writerObj.FrameRate = 25;
    open(writerObj);
    writeVideo(writerObj,frames2);
    close(writerObj);

end

% ret = RetargetMethods;
% path = 'F:\Tez\Thesis\Hollywood2-actions\Hollywood2\SubsetAVIS\Frames\actioncliptest00186.avi_Frames.mat';
% load(path);
% retFrames = ret.Cigdem(path,frames,[16*25 9*25]);


% ffmpeg conversion cmdline
% ffmpeg -i E:\Tez\Thesis\userStudy\aviwithgreen\actioncliptest00033.avi_Linear.avi.avi -c:v libx264 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -movflags faststart E:\Tez\Thesis\userStudy\aviwithgreen\actioncliptest00033_linear.mp4
