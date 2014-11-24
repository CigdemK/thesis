clear,clc;

folderName = 'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\';
fileName = 'actioncliptest00001.avi';
methodName = 'Cigdem';
shotStart = 48;
shotEnd = 95;

video=VideoReader([folderName fileName]);
frames = read(video);
frames = frames(:,:,:,shotStart:shotEnd);
ret = RetargetMethods;
tic;
switch methodName
    case 'Rubinstein2009'
        retFrames = ret.Rubinstein2009(frames,[round(video.Height*0.9) round(video.Width*0.9) ]);
    case 'Yan2013_Mine'
        retFrames = ret.Yan2013_Mine(frames,[round(video.Height*0.9) round(video.Width*0.9) ]);
    case 'Yan2013'
        retFrames = ret.Yan2013(frames,[round(video.Height*0.9) round(video.Width-20) ]);
    case 'Cigdem'
        retFrames = ret.Cigdem(frames,[round(video.Height*0.9) round(video.Width-20) ]);
    otherwise
        disp('No such retargeting algoritm!')
end
toc;
writerObj = VideoWriter([folderName fileName methodName '.avi']);
writerObj.FrameRate = video.FrameRate;
open(writerObj);tic;
writeVideo(writerObj,retFrames);toc;
close(writerObj);