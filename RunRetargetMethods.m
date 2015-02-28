clear,clc;

movieNames = {'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00001.avi'};
methodName = {'Wang'};
tic;
ret = RetargetMethods;
for i = 1:size(movieNames,2)
    
    video = VideoReader(movieNames{i});
    frames = read(video);
    
    for j = 1:size(methodName,2)
        
        currentMethod = methodName{j};
        switch currentMethod
            case 'Wang'
                retFrames = ret.Wang2011(frames,[16*20 9*20]);
            case 'Rubinstein2009'
                retFrames = ret.Rubinstein2009(frames,[video.Height 16*20]);save('tmp1.mat','retFrames');
                retFrames = ret.Rubinstein2009(retFrames,[9*20 16*20]);save('tmp2.mat','retFrames');
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
        
        writerObj = VideoWriter([movieNames{i} currentMethod '.avi']);
        writerObj.FrameRate = video.FrameRate;
        open(writerObj);
        writeVideo(writerObj,retFrames);
        close(writerObj);
        
    end
    toc;
end




% % save retargeted video
% writerObj = VideoWriter('F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00001.avi_mean.avi');
% writerObj.FrameRate = 23;
% open(writerObj);
% writeVideo(writerObj,newframes);
% close(writerObj);

% video = VideoReader('F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00001.aviCigdem.avi');
% retframes = read(video);
% 
% for i = 100
%     figure;
%     
%     imshow(newframes(:,:,:,i));
%     figure
%     imshow(retframes(:,:,:,i));
% end


