clc, clear,

% Read video
video = VideoReader('video000.avi');
nFrames = video.NumberOfFrames;
vidHeight = video.Height;
vidWidth = video.Width;

% Create movie structures for original, seam carved, linear scaled
mov(1:floor(nFrames)) = struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),'colormap',[]);
carvedMov(1:floor(nFrames)) = struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),'colormap', []);
scaledMov(1:floor(nFrames)) = struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),'colormap', []);

tic;
% Create original video
t = 1;
for k = 1 : nFrames
    mov(t).cdata = read(video, k);
    t = t+1;
end
toc;

% Create seam carved version of the original video
for k = 1 : (t-1)
    carvedMov(k).cdata = im2uint8(seamcarving(mov(k).cdata,200));
end
toc;

% Create linear scaled version of the original video
for k = 1 : (t-1)
    scaledMov(k).cdata = imresize(mov(k).cdata , [vidHeight vidWidth*3/4]);
end
toc;

% Play movie and save processed versions
movie2avi(mov, 'original.avi');
movie2avi(carvedMov, 'carved.avi');
movie2avi(scaledMov, 'scaled.avi');




