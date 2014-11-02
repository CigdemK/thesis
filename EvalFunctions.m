function F = EvalFunctions
    F.CalculateAUCscoreFrame = @CalculateAUCscoreFrame;
    F.CalculateAUCscoreVideo = @CalculateAUCscoreVideo;
end

function [ score , R ] = CalculateAUCscoreVideo( framesSaliency, framesFixation,algorithm )
    tic;
    nFrames = size(framesSaliency,2);

    R = zeros(12,2);
    for k = 1:nFrames
        
        [~,currentR] = CalculateAUCscoreFrame( framesSaliency{k}, framesFixation{k} );
        R = R + currentR;
        
        if mod(k, 100) == 0
            save(strcat('data\R_',algorithm,'_',num2str(k),'.mat'),'R');
            R = zeros(12,2);
        end
        
        toc;
    end
    score = 0;
%     R = R./nFrames;
%     score = trapz(flipdim(R(:,1),1),flipdim(R(:,2),1));
    
end

function [ score , R ] = CalculateAUCscoreFrame( salMap, eyeMap, shufMap, numRandom )

    if nargin < 3
        shufMap = true(size(eyeMap));
    end

    if nargin < 4
        numRandom = 50;
    end

    if isempty(shufMap) || max(max(shufMap)) == 0 % its empty or no fixation at all
        shufMap = true(size(eyeMap));
    end

    %%% Resize and normalize saliency map
    salMap = double(imresize(salMap,size(eyeMap),'nearest'));
    salMap = salMap - min(min(salMap));
    if max(max(salMap)) > 0
        salMap = salMap / max(max(salMap));
    end

    %%% Pick saliency value at each eye fixation along with [numrandom] random points
    [X Y] = find(eyeMap > 0);
    [XRest YRest] = find(shufMap > 0);
    localHum = nan(length(X),1);
    localRan = nan(length(X),numRandom);
    
    for k=1:length(X)
        localHum(k,1) = salMap(X(k),Y(k));
        for kk=1:numRandom
%             r = randi([1 length(XRest)],1);
            r = randi(length(XRest));
            localRan(k,kk) = salMap(XRest(r),YRest(r));
        end
    end

    %%% Calculate AUC score for each eyefix and randpoints
    R  = cell(1,size(localRan,2));
    steps = .1;
    asteps=0:steps:1;

    for ii = 1:size(localRan,2)
        R{ii} = auc_steps_helper(double(localHum),double(localRan(:, ii)),double(asteps));
    end

    x = 0;
    y = 0;
    for i = 1:size(R,2)
        x = x + R{i}(:,1);
        y = y + R{i}(:,2);
    end
    
    x=x/numRandom;
    y=y/numRandom;
    
    R = [x,y];
    score = trapz(flipdim(x,1),flipdim(y,1));
    
end



