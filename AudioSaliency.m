clc,clear;

filename = 'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\actioncliptest00001.avi';
load([filename '_videoSaliencyMap.mat']);
[vidHeight, vidWidth, nFrames] = size(videoSaliency);

a = mirfeatures(filename, 'Sampling' ,nFrames);
a = mirgetdata(a);

Util = UtilFunctions;
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


% canoncorr on seperate pixels
out = reshape(videoSaliency,[vidHeight*vidWidth,nFrames])';
[A, B, r, U ,V] = canoncorr(inp,out);
test = sum(B,2)';
for k = 1:nFrames
    out(k,:) = out(k,:).*test;
end
out = reshape(out, [vidHeight,vidWidth,nFrames]);
for k = 1:20:nFrames
    figure;
    imshow(mat2gray(out(:,:,k)))
end

% canoncorr on patches of 20x20
t = 1;
out = zeros(nFrames,ceil(vidWidth*vidHeight/400));
for i = 20:20:vidWidth
    for j = 20:20:vidHeight
        tmp = mean(mean(videoSaliency(j-19:j,i-19:i,:)));
        out(:,t) = tmp(:);
        t = t + 1;
    end
end
out(:,~any(out,1)) = [];
[A, B, r, U ,V] = canoncorr(inp,out);
test = sum(B,2)';
for k = 1:nFrames
    out(k,:) = out(k,:).*test;
end
out = reshape(out, [j/20,i/20,nFrames]);
for k = 1:20:nFrames
    figure;
    imshow(mat2gray(out(:,:,k)))
end











