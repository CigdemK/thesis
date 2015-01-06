clear,clc;

movieNames = {'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00001.avi'};
methodName = 'Cigdem';
tic;
ret = RetargetMethods;
for i = 1:size(movieNames,2)
    
    video = VideoReader(movieNames{i});
    frames = read(video);
    
    switch methodName
        case 'Rubinstein2009'
            retFrames = ret.Rubinstein2009(movieNames{i},[video.Height round(size(frames,2)-10)]);
        case 'Yan2013_Mine'
            retFrames = ret.Yan2013_Mine(frames,[round(video.Height) round(size(frames,2)-10)]);
        case 'Yan2013'
            retFrames = ret.Yan2013(frames,[video.Height round(size(frames,2)-50)]);
        case 'Linear'
            retFrames = ret.LinearScaling(frames,[304 516]);
        case 'Cigdem'
            retFrames = ret.Cigdem(movieNames{i},frames,[16*10 9*10]);
        otherwise
            disp('No such retargeting algoritm!')
    end
    toc;
end

writerObj = VideoWriter([movieNames{i} methodName '.avi']);
writerObj.FrameRate = video.FrameRate;
open(writerObj);tic;
writeVideo(writerObj,retFrames);
close(writerObj);






