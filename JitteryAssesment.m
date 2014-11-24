% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%  Yan, B., Yuan, B., and Yang, B. 2013.                 %
%  Effective Video Retargeting With Jittery Assessment.  %
%  IEEE Transactions on Multimedia, 11.                  %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %


% N       : nr of columns of the grid (x)
% M       : nr of rows of the grid (y)
% SRow    : Height of the row
% SColumn : Width of the column
% DRow    : Difference of rows of consecutive frames
% DColumn : Difference of columns of consecutive frames

clear, clc;

folderName = 'C:\Users\ckocb_000\Thesis\bin\others\';
fileName = 'rhinos.avi';
N = 20;
M = 10;

% Read the video
video=VideoReader([folderName fileName]);
nFrames = video.NumberOfFrames;
frames = read(video);

newSizeX = video.Width * 0.5; 
newSizeY = video.Height * 0.5;

SRow = zeros(M-1,nFrames);
SColumn = zeros(N-1,nFrames);
DRow = zeros(M-1,nFrames-1);
DColumn = zeros(N-1,nFrames-1);
jitterRow = zeros(nFrames-1,1);
jitterColumn = zeros(nFrames-1,1);

% Apply warping with a given mesh
WANG = WANGFunctions;
[~, meshX , meshY] = WANG.ScaleStretch(frames(:,:,:,1),[newSizeY, newSizeX], [M N]);

for column = 1:N-1
    SColumn(column,1) = meshX(1,column+1) - meshX(1,column);
end
for row = 1:M-1
    SRow(row,1) = meshY(row+1,1) - meshY(row,1);
end
tic;
for i = 2:10
    
    [~, meshX , meshY] = WANG.ScaleStretch(frames(:,:,:,i),[newSizeY, newSizeX], [M N]);
    
    for column = 1:N-1
        SColumn(column,i) = meshX(1,column+1) - meshX(1,column);
        DColumn(column,i-1) = SColumn(column,i) - SColumn(column,i-1);
    end
    for row = 1:M-1
        SRow(row,i) = meshY(row+1,1) - meshY(row,1);
        DRow(row,i-1) = SRow(row,i) - SRow(row,i-1);
    end

    toc;
end


for i = 1:20  
    jitterColumn(i) = sqrt(sum((DColumn(:,i).*1000).^2))/N;
    jitterRow(i) = sqrt(sum((DRow(:,i).*1000).^2))/M;
end

figure, title('Column Jitter'); hold on; plot(1:nFrames-1,jitterColumn);
figure, title('Row Jitter'); hold on; plot(1:nFrames-1,jitterRow);


 