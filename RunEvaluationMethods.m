clc,clear

videoPath = 'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00001.avi';
% methodName = {'Compress','PSNR','MSE','VIFP','VIF','VSNR','MSSIM','SSIM',...
%             'UQI','SNR','WSNR','NQM','IFC','DivisiveNormalization','Focus',...
%             'Blur','Sharpness','Brightness'};
methodName = {'Brightness','DivisiveNormalization'};

video=VideoReader(videoPath);
frames = read(video);
nFrames = size(frames,4);

Eval=EvaluationMetrics;
for i = 1:size(methodName,2)
    currentMethod = methodName{i};
    switch currentMethod
        case 'DivisiveNormalization'
            res = Eval.DivisiveNormalization(videoPath,frames);
        case 'Focus'
            res = Eval.Focus(videoPath,frames);
        case 'Blur'
            res = Eval.Blur(videoPath,frames);
        case 'Sharpness'
            res = Eval.Sharpness(videoPath,frames);
        case 'Brightness'
            res = Eval.Brightness(videoPath,frames);
        case 'MSE'
            res = Eval.MetrixMuxMSE(videoPath,frames);
        case 'PSNR'
            res = Eval.MetrixMuxPSNR(videoPath,frames);
        case 'VIFP'
            res = Eval.MetrixMuxVIFP(videoPath,frames);
        case 'VIF'
            res = Eval.MetrixMuxVIF(videoPath,frames);
        case 'VSNR'
            res = Eval.MetrixMuxVSNR(videoPath,frames);
        case 'MSSIM'
            res = Eval.MetrixMuxMSSIM(videoPath,frames);
        case 'SSIM'
            res = Eval.MetrixMuxSSIM(videoPath,frames);
        case 'UQI'
            res = Eval.MetrixMuxUQI(videoPath,frames);
        case 'SNR'
            res = Eval.MetrixMuxSNR(videoPath,frames);
        case 'WSNR'
            res = Eval.MetrixMuxWSNR(videoPath,frames);
        case 'NQM'
            res = Eval.MetrixMuxNQM(videoPath,frames);
        case 'IFC'
            res = Eval.MetrixMuxIFC(videoPath,frames);
        case 'Compress'
            res = Eval.CompressVideo(videoPath,frames);
        otherwise
            disp([methodName{i} ': No such evaluation method!']);
            res = [];
    end
    figure; plot(1:nFrames-1,res(1:end-1));
    xlabel('Frame Number','FontWeight','Bold');
    ylabel(methodName{i},'FontWeight','Bold');
    axis tight;
    
end
