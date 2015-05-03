clc,clear;
movieNames = {'E:\Hollywood2_part1\actioncliptest00001.avi',...
    'E:\Hollywood2_part1\actioncliptest00010.avi',...
    'E:\Hollywood2_part1\actioncliptest00020.avi',...
    'E:\Hollywood2_part1\actioncliptest00030.avi',...
    'E:\Hollywood2_part1\actioncliptest00040.avi',...
    'E:\Hollywood2_part1\actioncliptest00050.avi',...
    'E:\Hollywood2_part1\actioncliptest00060.avi',...
    'E:\Hollywood2_part1\actioncliptest00070.avi',...
    'E:\Hollywood2_part1\actioncliptest00080.avi',...
    'E:\Hollywood2_part1\actioncliptest00090.avi'};

I = ImageSaliency;
Util = UtilFunctions;
M = MotionFunctions;
nMovies = numel(movieNames);

Aall = zeros(28,9,nMovies);
Ball = zeros(9,9,nMovies);
rall = zeros(1,nMovies);
inpAll = cell(nMovies,1);
outAll = cell(nMovies,1);
tic;
for u = 1:3
    
    filename = movieNames{u};

%     videoSaliency = I.ActionsInTheEyeFixation(filename);
    video = VideoReader( filename );
    frames = read(video); 
    frames = frames(:,:,:,10:50); 
    [optY,optX] = M.MyOpticalFlow(frames);
    videoSaliency = (optY.^2 + optX.^2);
    [vidHeight, vidWidth, nFrames] = size(videoSaliency);

    a = mirfeatures(filename, 'Sampling' ,nFrames);
    a = mirgetdata(a);

    res = Util.RecursiveStructTraverse(a);

    nFeatures = numel(res);
    inp = zeros(nFrames,nFeatures);
    for i = 1:nFeatures
        res{i}(isnan(res{i})) = 0;
        if  numel(res{i}) > 9
            if numel(res{i}) > nFrames
                sampledRes{i} = downsample(res{i},round(size(res{i},2)/nFrames));
            elseif numel(res{i}) < nFrames
                sampledRes{i} = interp(res{i},round(nFrames/size(res{i},2)));
            else
                sampledRes{i} = res{i};
            end
            newSize = numel(sampledRes{i});
            if newSize < nFrames
                inp(:,i) = [sampledRes{i} repmat(sampledRes{i}(end), [1 abs(newSize-nFrames)])]';
            elseif newSize > nFrames
                inp(:,i) = sampledRes{i}(1:nFrames)';
            else
                inp(:,i) = sampledRes{i}';
            end
        end
    end

    % canoncorr on 3x3
    t = 1;
    out = zeros(nFrames,9);
    nHor = ceil(vidHeight/3);
    nVert = ceil(vidWidth/3);
    for i = 1:nVert:vidWidth
        for j = 1:nHor:vidHeight
            endx = min(i+nVert-1,vidWidth);
            endy = min(j+nHor-1,vidHeight);
            tmp = mean(mean(videoSaliency(j:endy,i:endx,:)));
            out(:,t) = tmp(:);
            t = t + 1;
        end
    end
%     out(:,~any(out,1)) = [];
%     inp(:,~any(inp,1)) = [];
    [A, B, r, U ,V] = canoncorr(inp,out);
    Aall(:,:,u) = A;
    Ball(:,:,u) = B;
    inpAll{u} = inp;
    outAll{u} = out;
    toc;
end



% audioFile = 'C:/Users/ckocb_000/Desktop/openSMILE-2.1.0/example-audio/opensmile.wav' ;
% [status,exeOutput] = system(['"C:/Users/ckocb_000/Desktop/openSMILE-2.1.0/SMILExtract_Release -C C:/Users/ckocb_000/Desktop/openSMILE-2.1.0/config/MFCC12_E_D_A.conf -I C:/Users/ckocb_000/Desktop/openSMILE-2.1.0/example-audio/opensmile.wav"']);
% outputFile = 'C:/Users/ckocb_000/Desktop/openSMILE-2.1.0/outputs/output_file.csv';
% mfcfile = fopen( outputFile, 'r', 'b' );
% 
% nSamples = fread( mfcfile, 1, 'int32' );
% sampPeriod = fread( mfcfile, 1, 'int32' )*1E-7;
% sampSize = 0.25*fread( mfcfile, 1, 'int16' );
% parmKind = fread( mfcfile, 1, 'int16' );
% 
% features = fread( mfcfile, [ sampSize, nSamples ], 'float' ).';
% 
% fclose( mfcfile );


for i = 1:1
corr(Ball(:,1,i),mean(outAll{i})')
end
% save('audioPaper/groundTruth.mat','Aall','Ball','inpAll','outAll')







