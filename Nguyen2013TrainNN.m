%% INITIALIZATION

clc, clear;

% Read fixation maps
pathStatic  = 'C:\Users\ckocb_000\Thesis\BACKUP\CAMO\Static-Fixation\';
pathDynamic = 'C:\Users\ckocb_000\Thesis\BACKUP\CAMO\Dynamic-Fixation\';
pathKeyFrames = 'C:\Users\ckocb_000\Thesis\BACKUP\CAMO\AVISelectedVideos_KeyFrames\';
pathPrevKeyFrames = 'C:\Users\ckocb_000\Thesis\BACKUP\CAMO\AVISelectedVideos_PreviousKeyFrames\';

fileName  = dir(fullfile(pathStatic, '*.jpg')); % all dir paths includes files with the same name, lets use pathStatic
fileNamePrev  = dir(fullfile(pathPrevKeyFrames, '*.jpg')); % all dir paths includes files with the same name, lets use pathStatic

C = cell(length(fileName), 2);
for k = 1:length(fileName)
    currentFileName = fileName(k).name;
    currentFileNamePrev = fileNamePrev(k).name;
    C{k,1} = imread(strcat(pathStatic,currentFileName));
    C{k,2} = imread(strcat(pathDynamic,currentFileName));
    C{k,3} = imread(strcat(pathKeyFrames,currentFileName));
    C{k,4} = imread(strcat(pathPrevKeyFrames,currentFileNamePrev));
end

% Init variables
[ height width ] = size( C{ 1 , 1 } );
numberOfPatches = length( C ) * floor( height / 90 ) * floor( width / 90 );
p  = zeros(90,90,numberOfPatches);
xy = zeros(numberOfPatches,2);
pi = zeros(90,90,numberOfPatches);
pv = zeros(90,90,numberOfPatches);
CM = zeros(3,3,numberOfPatches);

% set optical flow parameters (see Coarse2FineTwoFrames.m for the definition of the parameters)
alpha = 0.012;
ratio = 0.75;
minWidth = 20;
nOuterFPIterations = 7;
nInnerFPIterations = 1;
nSORIterations = 30;

para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];

%% PREPARE PATCHES

i=1;
for n = 1 : length( C )
    tic;
    dynamicFixation = C{n,2}; %take the dynamic fixation map as the ground truth
    currentImage = C{n,3};
    previousImage = C{n,4};
    
%     Calculate dynamic and static saliency 
    [vx,vy,warpI2] = Coarse2FineTwoFrames(previousImage,currentImage,para);
    currentVideo = sqrt(vx.^2+vy.^2);
    currentImage = saliency(currentImage); 
    
    % Divide into patches
    for k = 1 : floor( height / 90 )
        for t = 1 : floor( width / 90 )
            
            p( :,:,i )  = dynamicFixation( (k*90-89) : (k*90) , (t*90-89) : (t*90) );
            pi( :,:,i ) = currentImage( (k*90-89) : (k*90) , (t*90-89) : (t*90) );
            pv( :,:,i ) = currentVideo( (k*90-89) : (k*90) , (t*90-89) : (t*90));
            xy(i,1) = k;
            xy(i,2) = t;
           
%             meanX = mean(mean(vx((k*90-89) : (k*90) , (t*90-89) : (t*90))));
%             meanY = mean(mean(vy((k*90-89) : (k*90) , (t*90-89) : (t*90))));
%             vxTemp = mean(mean([(vx((k*90-89) : (k*90) , (t*90-89) : (t*90))<meanX*0.9) (vx((k*90-89) : (k*9) , (t*9-8) : (t*9))>meanX*(-0.9))]));
%             vyTemp = mean(mean([(vy((k*90-89) : (k*90) , (t*90-89) : (t*90))<meanY*0.9) (vy((k*90-89) : (k*9) , (t*9-8) : (t*9))>meanY*(-0.9))]));

            i = i + 1;
        end
    end
    toc;
end


%% TRAIN NEURAL NETWORK

one   =CM(1,1,:);
two   =CM(1,2,:);
three =CM(1,3,:);
four  =CM(2,1,:);
five  =CM(2,2,:);
six   =CM(2,3,:);
seven =CM(3,1,:);
eight =CM(3,2,:);
nine  =CM(3,3,:);
inputs = [ one(:)   , two(:)   , three(:) , ...
           four(:)  , five(:)  , six(:)   , ...
           seven(:) , eight(:) , nine(:)  , ...
           xy(:,1)  , xy(:,2)]'; % input vector (6-dimensional pattern)
       
wi = repmat(0.5,[81, numberOfPatches]) ; % corresponding target output vector
wv = repmat(0.5,[81, numberOfPatches]) ; % corresponding target output vector

neti = network( ...
    1, ... % numInputs, number of inputs,
    3, ... % numLayers, number of layers
    [1;1; 1], ... % biasConnect, numLayers-by-1 Boolean vector,
    [1;0; 0], ... % inputConnect, numLayers-by-numInputs Boolean matrix,
    [0 0 0; 1 0 0 ; 0 1 0], ... % layerConnect, numLayers-by-numLayers Boolean matrix
    [0 0 1] ... % outputConnect, 1-by-numLayers Boolean vector
    );

netv = network( ...
    1, ... % numInputs, number of inputs,
    3, ... % numLayers, number of layers
    [1;1; 1], ... % biasConnect, numLayers-by-1 Boolean vector,
    [1;0; 0], ... % inputConnect, numLayers-by-numInputs Boolean matrix,
    [0 0 0; 1 0 0 ; 0 1 0], ... % layerConnect, numLayers-by-numLayers Boolean matrix
    [0 0 1] ... % outputConnect, 1-by-numLayers Boolean vector
    );

% number of hidden layer neurons
neti.layers{1}.size = 5;
neti.layers{2}.size = 1;
netv.layers{1}.size = 5;
netv.layers{2}.size = 1;

% hidden layer transfer function
neti.layers{1}.transferFcn = 'tansig';
neti.layers{2}.transferFcn = 'purelin';
netv.layers{1}.transferFcn = 'tansig';
netv.layers{2}.transferFcn = 'purelin';


% network training
neti.trainFcn = 'trainlm';
neti.performFcn = 'mse';
netv.trainFcn = 'trainlm';
netv.performFcn = 'mse';

fi = configure(neti,inputs(:,1),wi(:,1));
fv = configure(netv,inputs(:,1),wv(:,1));
fi = train(neti,inputs(:,1),wi(:,1));
fv = train(netv,inputs(:,1),wv(:,1));

for j = 1 : 4200
    
    t=0; 
    lambda = 0.1;
    ro = 1.5;
    tic;
    while(t<30)
        t=t+1;
        fiOld = fi;
        fi = configure(fi,inputs(:,j),wi(:,t+1));
        fv = configure(fv,inputs(:,j),wv(:,t));
        
        wi_90 = reshape(repmat(wi(:,t),10),[90 90]); % ameleligin alasi
        wv_90 = reshape(repmat(wv(:,t),10),[90 90]); % ameleligin alasi
        wi_last = reshape(repmat(lambda * ( fiOld(inputs(:,j)) + fv(inputs(:,j)) - wv(: ,t) ) , 10),[90 90]);
        wv_last = reshape(repmat(lambda * ( fiOld(inputs(:,j)) + fv(inputs(:,j)) - wi(: ,t) ) , 10),[90 90]);
        
        witmp = ( pi( : , :  , j )' * pi( : , : , j ) + lambda' * eye( 90 ) ) \  ...
                   ( pi( : , : , j ) * p( : , : , j )' - pi( : , : , j) * pv( : , : , j )' * (wi_90) ...
                   + wi_last );
        wvtmp = ( pv( : , : , j )' * pv( : , : , j ) + lambda' * eye( 90 ) ) \  ...
                   ( pv( : , : , j ) * p( : , : , j )' - pv( : , : , j) * pi( : , : , j )' * (wv_90) ...
                   + wv_last);
               
        i=1;       
        for k =1:10:90
            for s = 1:10:90
                wi(i,t+1) = mean(mean(witmp(k:k+9,s:s+9)));
                wv(i,t+1) = mean(mean(wvtmp(k:k+9,s:s+9)));
                i = i+1;
            end
        end
        
        fi = train(fi,inputs(:,j),wi(:,t+1));
        fv = train(fv,inputs(:,j),wv(:,t+1));

        wi(:,t+1) = fi(inputs(:,j));
        wv(:,t+1) = fv(inputs(:,j));

        lambda = ro * lambda;
    end
    toc;
end


%% CREATE A CROPPED MOVIE WITH THE TRAINED NETWORK
clear;clc;
load('f_w.mat');
FolderName = '../CAMO_Videos/Tilt/';
FileName = 'Tilt0015.avi';


RATIO = [16 9]; % aspect ratio of smartphones
CROP =10; %the crop value will be multiplied with aspect ration to get crop window

CROPPINGX = CROP * RATIO(1,1);
CROPPINGY = CROP * RATIO(1,2);

video=VideoReader([folderName fileName]);
frames = read(video);

% set optical flow parameters (see Coarse2FineTwoFrames.m for the definition of the parameters)
alpha = 0.012;
ratio = 0.75;
minWidth = 20;
nOuterFPIterations = 7;
nInnerFPIterations = 1;
nSORIterations = 30;

para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];


% Create cropped movie with cropX cropY values
tic;
energyMap = zeros(video.Height,video.Width,video.NumberOfFrames);
opticalFlowMap = zeros(video.Height,video.Width,video.NumberOfFrames);
saliencyMap = zeros(video.Height,video.Width,video.NumberOfFrames);
grayFrames = zeros(video.Height,video.Width,video.NumberOfFrames);


for n = 1 : video.NumberOfFrames-1
    im1 = frames(:,:,:,n);
    im2 = frames(:,:,:,n+1);
    [vx,vy,warpI2] = Coarse2FineTwoFrames(im1,im2,para);
    toc;
    dynamicSaliency = sqrt(vx.^2+vy.^2);
    staticSaliency = saliency(im1); 
    
    opticalFlowMap(:,:,n) = dynamicSaliency;
    saliencyMap(:,:,n) = staticSaliency;
    toc;
    try
        cd('../CameraMotionCode');
        [h mappedImg2] = homography( im1,im2);
        cd('../bin');
    catch err
        h = zeros(3);   
    end
    toc;
    currentPatch = [];
    % Divide into patches
    i=1;
    nrPatches = floor( video.Height / 90 )*floor( video.Width / 90 );
    xyTest = zeros(nrPatches,2);
    CM = zeros(3,3,nrPatches);
    currentFrame=zeros(video.Height,video.Width);
    for k = 1 : floor( video.Height / 90 )
        for t = 1 : floor( video.Width / 90 )
            xyTest(i,1) = k;
            xyTest(i,2) = t;
            CM(:,:,i) = h;
            
            one   =CM(1,1,i);
            two   =CM(1,2,i);
            three =CM(1,3,i);
            four  =CM(2,1,i);
            five  =CM(2,2,i);
            six   =CM(2,3,i);
            seven =CM(3,1,i);
            eight =CM(3,2,i);
            nine  =CM(3,3,i);
            inputs = [ one(:)   , two(:)   , three(:) , ...
                       four(:)  , five(:)  , six(:)   , ...
                       seven(:) , eight(:) , nine(:)  , ...
                       xyTest(i,1)  , xyTest(i,2)]'; % input vector (6-dimensional pattern)
                          
            wi = fi(inputs);
            wv = fv(inputs);
            wi = reshape(repmat(wi,10),[90 90]);
            wv = reshape(repmat(wv,10),[90 90]);
            
            videoSaliency = wi*staticSaliency((k*90-89) : (k*90) , (t*90-89) : (t*90))  + ...
                wv*dynamicSaliency((k*90-89) : (k*90) , (t*90-89) : (t*90) );
            currentFrame(k*90-89:k*90,t*90-89:t*90) = videoSaliency(:,:);
                        
            i = i + 1;
        end
    end
    
    energyMap(:,:,n) = currentFrame;
    grayFrames(:,:,n) = mat2gray(energyMap(:,:,n));
    toc;  
end

% save(strcat(FolderName,FileName,'.mat'),'grayFrames','energyMap','saliencyMap','dynamicSaliency');

