clc,clear;

warning off MATLAB:MKDIR:DirectoryExists

F = dir('E:\Tez\Thesis\userStudy\originalAvi\*.avi');

methodName = {'mySSIM'};

tic;
Eval=EvaluationMetrics;
for j = 1:length(F)

    D = struct('methods',{},'results',[]);
    D(1).methods = {0};
    D(1).results = {0};

    videoPath = ['E:\Tez\Thesis\userStudy\originalAvi\' F(j).name];
    
    res = regexp(F(j).name,'_','split');
    originalPath = ['E:\Tez\Thesis\userStudy\originalAvi\linear\' res{1} '_Linear.avi'];
    
    video = VideoReader(videoPath);
    frames = read(video);
    
    video = VideoReader(originalPath);
    originalFrames = read(video);
    
    for i = 1:size(methodName,2)
        
        currentMethod = methodName{i};
        
        switch currentMethod
            case 'JerkinessMatlab'
                res = Eval.JerkinessMatlab(videoPath,frames);
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
            case 'mySSIM'
                res = Eval.SSIM(originalFrames,frames);
            otherwise
                disp([methodName{i} ': No such evaluation method!']);
                res = [];
        end
        D.methods{end+1} = currentMethod;
        D.results{end+1} = res;
    end
    save(['E:\Tez\Thesis\userStudy\metricResults\' F(j).name '_ssim.mat'],'D');
    toc;
end



% methodName = {'Compress','PSNR','MSE','VIFP','VIF','VSNR','MSSIM','SSIM',...
%             'UQI','SNR','WSNR','NQM','IFC','DivisiveNormalization','Focus',...
%             'Blur','Sharpness','Brightness'};
% {'Focus','Blur','Sharpness','Brightness','UQI','Compress','MSE','PSNR'}



