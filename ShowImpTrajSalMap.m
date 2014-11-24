clc,clear;
folderName = 'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\';
fileName = 'actioncliptest00340.avi'; %_ImprovedTrajectoriesSaliencyMap.mat

video=VideoReader([folderName fileName]);
frames = read(video);

load('F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00340\ImprovedTrajectoriesSaliencyMap.mat');

for i = 1:20:video.NumberOfFrames
    figure;
    subplot(2,1,1);
    imshow(frames(:,:,:,i));
    subplot(2,1,2);
    imshow(saliencyMap(:,:,i));
end