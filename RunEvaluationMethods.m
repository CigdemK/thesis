clc,clear;

warning off MATLAB:MKDIR:DirectoryExists
load('D.mat');
for j = 91:870
    currentRes = struct(); 
    currentRes.path = ['F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest000' num2str(j+9) '.avi'];
    currentRes.methods = cell(1,1);
    currentRes.results = zeros(1,1);
    D(j) = currentRes;
end

methodName = {'UQI'};
%100,70,130,62,10,13,34,29,18
tic;
Eval=EvaluationMetrics;
for j = 91:870
    
    videoPath = ['F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00' num2str(j+9) '.avi'];
    
    video = VideoReader(videoPath);
    frames = read(video);
    
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
            otherwise
                disp([methodName{i} ': No such evaluation method!']);
                res = [];
        end

        D(j).methods{end+1} = currentMethod;
        if length(res)>1
            D(j).results(end+1) = mean(res(1:end-1));
        else
            D(j).results(end+1) = res;
        end
        toc;
    end
    
end

% Dnew = D;
% load('D.mat');
% for i = 1: size(D,2)
%     D(i).methods = [D(i).methods,Dnew(i).methods(2:end)];
%     D(i).results = [D(i).results,Dnew(i).results(2:end)];
% end
save('Dnew.mat','D');

max = 0.885783403663150;
ind = 0;
for i = 91:112
    if D(i).results(2)>max
        max = D(i).results(2);
        ind = i;
    end
end
% % delete info from D
% clc,clear;
% load('D.mat');
% for i = 1: size(D,2)
%     D(i).methods = D(i).methods(2:end);
%     D(i).results = D(i).results(2:end);
% end
% save('D.mat','D');


% methodName = {'Compress','PSNR','MSE','VIFP','VIF','VSNR','MSSIM','SSIM',...
%             'UQI','SNR','WSNR','NQM','IFC','DivisiveNormalization','Focus',...
%             'Blur','Sharpness','Brightness'};
% {'Focus','Blur','Sharpness','Brightness','UQI','Compress','MSE','PSNR'}