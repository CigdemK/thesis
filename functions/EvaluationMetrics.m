function F = EvaluationMetrics 
    F.CompressVideo = @CompressVideo;
    F.DivisiveNormalization = @DivisiveNormalization;
    F.MetrixMuxMSE = @MetrixMuxMSE;
    F.MetrixMuxPSNR = @MetrixMuxPSNR;
    F.MetrixMuxVIFP = @MetrixMuxVIFP;
    F.MetrixMuxVIF = @MetrixMuxVIF;
    F.MetrixMuxVSNR = @MetrixMuxVSNR;
    F.MetrixMuxMSSIM = @MetrixMuxMSSIM;
    F.MetrixMuxSSIM = @MetrixMuxSSIM;
    F.MetrixMuxUQI = @MetrixMuxUQI;
    F.MetrixMuxSNR = @MetrixMuxSNR;
    F.MetrixMuxWSNR = @MetrixMuxWSNR;
    F.MetrixMuxNQM = @MetrixMuxNQM;
    F.MetrixMuxIFC = @MetrixMuxIFC;

end

function D = MetrixMuxPSNR(videoPath,frames) % peak signal-to-noise ratio 

    if nargin < 2
        video=VideoReader(videoPath);
        frames = read(video);
    end
    
    nFrames = size(frames,4);
    
    D = zeros(nFrames,1);
    for i = 1:nFrames-1
        D(i) = metrix_mux(frames(:,:,:,i),frames(:,:,:,i+1),'PSNR');
    end
    
end

function D = MetrixMuxMSE(videoPath,frames) % mean-squared error 

    if nargin < 2
        video=VideoReader(videoPath);
        frames = read(video);
    end
    
    nFrames = size(frames,4);
    
    D = zeros(nFrames,1);
    for i = 1:nFrames-1
        D(i) = metrix_mux(frames(:,:,:,i),frames(:,:,:,i+1),'MSE');
    end
    
end

function D = MetrixMuxVIFP(videoPath,frames) % pixel-based VIF  

    if nargin < 2
        video=VideoReader(videoPath);
        frames = read(video);
    end
    
    nFrames = size(frames,4);
    
    D = zeros(nFrames,1);
    for i = 1:nFrames-1
        D(i) = metrix_mux(frames(:,:,:,i),frames(:,:,:,i+1),'VIFP');
    end
    
end

function D = MetrixMuxVIF(videoPath,frames) % visual information fidelity

    if nargin < 2
        video=VideoReader(videoPath);
        frames = read(video);
    end
    
    nFrames = size(frames,4);
    
    D = zeros(nFrames,1);
    for i = 1:nFrames-1
        D(i) = metrix_mux(frames(:,:,:,i),frames(:,:,:,i+1),'VIF');
    end
    
end

function D = MetrixMuxVSNR(videoPath,frames) % visual signal-to-noise ratio

    if nargin < 2
        video=VideoReader(videoPath);
        frames = read(video);
    end
    
    nFrames = size(frames,4);
    
    D = zeros(nFrames,1);
    for i = 1:nFrames-1
        D(i) = metrix_mux(frames(:,:,:,i),frames(:,:,:,i+1),'VSNR');
    end
    
end

function D = MetrixMuxMSSIM(videoPath,frames) % multiscale SSIM index 

    if nargin < 2
        video=VideoReader(videoPath);
        frames = read(video);
    end
    
    nFrames = size(frames,4);
    
    D = zeros(nFrames,1);
    for i = 1:nFrames-1
        D(i) = metrix_mux(frames(:,:,:,i),frames(:,:,:,i+1),'MSSIM');
    end
    
end

function D = MetrixMuxSSIM(videoPath,frames) % structural similarity index

    if nargin < 2
        video=VideoReader(videoPath);
        frames = read(video);
    end
    
    nFrames = size(frames,4);
    
    D = zeros(nFrames,1);
    for i = 1:nFrames-1
        D(i) = metrix_mux(frames(:,:,:,i),frames(:,:,:,i+1),'SSIM');
    end
    
end

function D = MetrixMuxUQI(videoPath,frames) % universal quality index

    if nargin < 2
        video=VideoReader(videoPath);
        frames = read(video);
    end
    
    nFrames = size(frames,4);
    
    D = zeros(nFrames,1);
    for i = 1:nFrames-1
        D(i) = metrix_mux(frames(:,:,:,i),frames(:,:,:,i+1),'UQI');
    end
    
end

function D = MetrixMuxIFC(videoPath,frames) % image fidelity criterion

    if nargin < 2
        video=VideoReader(videoPath);
        frames = read(video);
    end
    
    nFrames = size(frames,4);
    
    D = zeros(nFrames,1);
    for i = 1:nFrames-1
        D(i) = metrix_mux(frames(:,:,:,i),frames(:,:,:,i+1),'IFC');
    end
    
end

function D = MetrixMuxNQM(videoPath,frames) % noise quality measure

    if nargin < 2
        video=VideoReader(videoPath);
        frames = read(video);
    end
    
    nFrames = size(frames,4);
    
    D = zeros(nFrames,1);
    for i = 1:nFrames-1
        D(i) = metrix_mux(frames(:,:,:,i),frames(:,:,:,i+1),'NQM');
    end
    
end

function D = MetrixMuxWSNR(videoPath,frames) % weighted signal-to-noise ratio

    if nargin < 2
        video=VideoReader(videoPath);
        frames = read(video);
    end
    
    nFrames = size(frames,4);
    
    D = zeros(nFrames,1);
    for i = 1:nFrames-1
        D(i) = metrix_mux(frames(:,:,:,i),frames(:,:,:,i+1),'WSNR');
    end
    
end

function D = MetrixMuxSNR(videoPath,frames) % signal-to-noise ratio

    if nargin < 2
        video=VideoReader(videoPath);
        frames = read(video);
    end
    
    nFrames = size(frames,4);
    
    D = zeros(nFrames,1);
    for i = 1:nFrames-1
        D(i) = metrix_mux(frames(:,:,:,i),frames(:,:,:,i+1),'SNR');
    end
    
end

function D = DivisiveNormalization(videoPath,frames)

    if nargin < 2
        video=VideoReader(videoPath);
        frames = read(video);
    end
   
    nFrames = size(frames,4);
    
    D = zeros(nFrames,1);
    for i = 1:nFrames-1
        [d,~,~] = div_norm_metric(frames(:,:,:,i),frames(:,:,:,i+1));
        D(i) = d;
    end
    
end

function compressRatio = CompressVideo(videoPath,frames)

    if nargin < 2
        video=VideoReader(videoPath);
        frames = read(video);
    end

    originalFileProp = dir(videoPath);
    originalSize = originalFileProp.bytes;
    
    regexres = regexp(videoPath,'.avi','split');
    newVideoPath = [regexres{1} '_mpeg.mp4'];
    
    writerObj = VideoWriter(newVideoPath,'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,frames);
    close(writerObj);
    
    retFileProp = dir(newVideoPath);
    retSize = retFileProp.bytes;
    delete(newVideoPath);
    
    compressRatio = ( 100 * retSize ) / originalSize;
    
end