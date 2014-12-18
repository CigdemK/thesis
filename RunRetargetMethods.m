clear,clc;

folderName = 'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\';
fileName = 'actioncliptest00001.avi';
methodName = 'Cigdem';

video=VideoReader([folderName fileName]);
frames = read(video);

ret = RetargetMethods;
switch methodName
    case 'Rubinstein2009'
        for i = 1:5
%             load([folderName fileName 'Rubinstein2009.mat'],'retFrames');
%             frames = retFrames;
            retFrames = ret.Rubinstein2009(retFrames,[video.Height round(size(frames,2)-10)]);
            save([folderName fileName 'Rubinstein2009.mat'],'retFrames');
        end
    case 'Yan2013_Mine'
        load([folderName fileName '_Yan2013_Mine.mat'],'retFrames');
        frames = retFrames;
        retFrames = ret.Yan2013_Mine(frames,[round(video.Height) round(size(frames,2)-10)]);
        save([folderName fileName '_Yan2013_Mine.mat'],'retFrames');
    case 'Yan2013'
        retFrames = ret.Yan2013(frames,[video.Height round(size(frames,2)-50)]);
        save([folderName fileName '_Yan2013.mat'],'retFrames');
    case 'Linear'
        retFrames = ret.LinearScaling(frames,[200 355]);
        save([folderName fileName 'LinearScaling.mat'],'retFrames');
    case 'Cigdem'
        retFrames = ret.Cigdem(frames,[9*25 16*25]);
    otherwise
        disp('No such retargeting algoritm!')
end

writerObj = VideoWriter([folderName fileName methodName '.avi']);
writerObj.FrameRate = video.FrameRate;
open(writerObj);tic;
writeVideo(writerObj,retFrames);
close(writerObj);