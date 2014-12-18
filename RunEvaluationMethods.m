clc,clear

warning off MATLAB:MKDIR:DirectoryExists

videoPath = 'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00001.aviLinear.avi';
methodName = {'Shakiness'};
% methodName = {'Focus','Blur','Sharpness','Brightness','Compress','PSNR','MSE', ...
%             'VIFP','VSNR','MSSIM','SSIM','UQI','SNR','WSNR'};
% methodName = {'Compress','PSNR','MSE','VIFP','VIF','VSNR','MSSIM','SSIM',...
%             'UQI','SNR','WSNR','NQM','IFC','DivisiveNormalization','Focus',...
%             'Blur','Sharpness','Brightness'};

video=VideoReader(videoPath);
frames = read(video);

frames = frames(:,:,:,1:10);


nFrames = size(frames,4);
D=cell(1,1);
Eval=EvaluationMetrics;
for i = 1:size(methodName,2)
    currentMethod = methodName{i};
    switch currentMethod
        case 'Shakiness'
            res = Eval.Shakiness(videoPath,frames);
        case 'DivisiveNormalization'
            res = Eval.DivisiveNormalization(videoPath,frames);
        case 'Focus'
            res = Eval.Focus(videoPath,frames);
        case 'Blur'
            res = Eval.Blur(videoPath,frames);
        case 'MDE'
            res = Eval.MDE(videoPath,frames);
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
    
    D{end+1} = {methodName{i},res}
    
%     figure; plot(1:nFrames-1,res(1:end-1));
%     xlabel('Frame Number','FontWeight','Bold');
%     ylabel(methodName{i},'FontWeight','Bold');
%     axis tight;
    
end

