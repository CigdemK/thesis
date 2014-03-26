% % clc, clear;

% TRAIN NEURAL NETWORK

load('C:\Users\ckocb_000\Thesis\bin\H.mat');

one   =H(1,1,:);
two   =H(1,2,:);
three =H(1,3,:);
four  =H(2,1,:);
five  =H(2,2,:);
six   =H(2,3,:);
seven =H(3,1,:);
eight =H(3,2,:);

% one   =H(1,1,50);
% two   =H(1,2,50);
% three =H(1,3,50);
% four  =H(2,1,50);
% five  =H(2,2,50);
% six   =H(2,3,50);
% seven =H(3,1,50);
% eight =H(3,2,50);

outputs = zeros(6,length(inputs));
for i = 1:length(inputs)
    class = M(i);
    if      i < 26  class = 1;
    elseif  i < 52  class = 2;
    elseif  i < 78  class = 3;
    elseif  i < 114 class = 4;
    elseif  i < 140 class = 5;
    else            class = 6;
    end

    outputs(class,i) = 1;
end
       
inputs = [one(:)   , two(:)   , three(:) , ...
           four(:)  , five(:)  , six(:)   , ...
           seven(:) , eight(:) ]';
[inputs,PS] = mapminmax(inputs,0,+1);

neti = network( ...
    1, ... % numInputs, number of inputs,
    2, ... % numLayers, number of layers
    [1;1], ... % biasConnect, numLayers-by-1 Boolean vector,
    [1;0], ... % inputConnect, numLayers-by-numInputs Boolean matrix,
    [0 0; 1 0 ], ... % layerConnect, numLayers-by-numLayers Boolean matrix
    [0 1] ... % outputConnect, 1-by-numLayers Boolean vector
    );

% number of hidden layer neurons
neti.layers{1}.size = 1;
% neti.layers{2}.size = 1;

% hidden layer transfer function
neti.layers{1}.transferFcn = 'purelin';
% neti.layers{2}.transferFcn = 'transig';

% network training
neti.trainFcn = 'trainlm';
neti.performFcn = 'mse';

neti.trainParam.epochs = 100000;

neti = configure(neti,inputs,outputs);
neti = train(neti,inputs,outputs);


%% READ DATA FROM MOVIES

FolderName = '../CAMO_Videos/Tilt/';
FileName = 'Tilt0014.avi';
UtilFunctions = Functions;
UtilFunctions.ReadData(FolderName,FileName);
load('Data.mat');
mov = UtilFunctions.NewMovie(nFrames , vidHeight   ,vidWidth);
mov = UtilFunctions.ReadMovie(mov , video , nFrames );


H= [];
for n = 1 : nFrames-1
    tic;
    im1 = mov(n).cdata;
    im2 = mov(n+1).cdata;
    try
        cd('../CameraMotionCode');
        [h mappedImg2] = homography( im1 , im2 );
        cd('../bin');
    catch err
        h = zeros(3);   
    end
    H = cat(3,H,h);
    toc;
    n
end

save(strcat(FolderName,FileName,'.mat'),'H');

%% CREATE INPUT DATA

FolderName = '../CAMO_Videos/Tilt/';
FileName = 'Tilt0014.avi';
load(strcat(FolderName,FileName,'.mat'));

% zoomInput = [];
for i = 1:5:length(H)-4
    currentInput = [];
    for j = 1:4
        currentInput = [currentInput ; H(1,1,i+j) ; H(1,2,i+j) ; ...
                        H(1,3,i+j); H(2,1,i+j) ; H(2,2,i+j) ; ...
                        H(2,3,i+j); H(3,1,i+j) ; H(3,2,i+j)];
    end
    tiltInput = [ tiltInput , currentInput ]; 
end

% inputs = [dollyInput(:,1:26) , panInput(:,1:26) , pedestalInput(:,1:26) ...
%           tiltInput(:,1:26) , truckingInput(:,1:26) , zoomInput(:,1:26) ];

save('inputs.mat','inputs', 'dollyInput' , 'panInput' , ...
    'pedestalInput' ,'tiltInput' , 'truckingInput' , 'zoomInput');
