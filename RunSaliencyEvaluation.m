clear,clc;

movieNames = {'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00010.avi',...
            'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00004.avi'};
%             'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00001.avi',...
%             'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00031.avi',...
%             'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00084.avi'
        
Eval = EvalFunctions;
RImp = zeros(12,2);
RJudd = zeros(12,2);
RItti = zeros(12,2);
tic;
for i = 1:size(movieNames,2)

    AITE = load([movieNames{i} '_AITE_Saliency.mat']);
    Imp = load([movieNames{i} '_ImpTrj_Saliency.mat']);
    Judd = load([movieNames{i} '_staticSaliencyMap.mat']);
    Itti = load([movieNames{i} '_Itti_Saliency.mat']);

    nFrames = size(Imp.saliencyMap,3);
    
    ind = [];
    for k = 1:nFrames
        if sum(sum(Imp.saliencyMap(:,:,k)))~=0
            ind = [ind;k];
        end
    end
    
%     [~, R] = Eval.CalculateAUCscoreVideo(Imp.saliencyMap(:,:,ind),AITE.saliencyMap(:,:,ind));
%     RImp = RImp + R;
%     toc;
%     [~, R] = Eval.CalculateAUCscoreVideo(Judd.saliencyMap(:,:,ind),AITE.saliencyMap(:,:,ind));
%     RJudd = RJudd + R;
%     toc;
    [~, R] = Eval.CalculateAUCscoreVideo(Itti.saliencyMap(:,:,ind),AITE.saliencyMap(:,:,ind))
    RItti = RItti + R;
    toc;
    
end

% RImp  = RImp./size(movieNames,2);
% RJudd = RJudd./size(movieNames,2);
RItti = RItti./size(movieNames,2);

% scoreImp  = trapz(flipdim(RImp(:,1),1),flipdim(RImp(:,2),1))
% scoreJudd = trapz(flipdim(RJudd(:,1),1),flipdim(RJudd(:,2),1))
scoreItti = trapz(flipdim(RItti(:,1),1),flipdim(RItti(:,2),1))







