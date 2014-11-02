% % % clc, clear;
% 
% % TRAIN NEURAL NETWORK
% 
% load('C:\Users\ckocb_000\Thesis\bin\H.mat');
% 
% one   =H(1,1,:);
% two   =H(1,2,:);
% three =H(1,3,:);
% four  =H(2,1,:);
% five  =H(2,2,:);
% six   =H(2,3,:);
% seven =H(3,1,:);
% eight =H(3,2,:);
% 
% % one   =H(1,1,50);
% % two   =H(1,2,50);
% % three =H(1,3,50);
% % four  =H(2,1,50);
% % five  =H(2,2,50);
% % six   =H(2,3,50);
% % seven =H(3,1,50);
% % eight =H(3,2,50);
% 
% outputs = zeros(6,length(inputs));
% for i = 1:length(inputs)
%     class = M(i);
%     if      i < 26  class = 1;
%     elseif  i < 52  class = 2;
%     elseif  i < 78  class = 3;
%     elseif  i < 114 class = 4;
%     elseif  i < 140 class = 5;
%     else            class = 6;
%     end
% 
%     outputs(class,i) = 1;
% end
% 
%        
% inputs = [one(:)   , two(:)   , three(:) , ...
%            four(:)  , five(:)  , six(:)   , ...
%            seven(:) , eight(:) ]';
% [inputs,PS] = mapminmax(inputs,0,+1);
% 
% neti = network( ...
%     1, ... % numInputs, number of inputs,
%     2, ... % numLayers, number of layers
%     [1;1], ... % biasConnect, numLayers-by-1 Boolean vector,
%     [1;0], ... % inputConnect, numLayers-by-numInputs Boolean matrix,
%     [0 0; 1 0 ], ... % layerConnect, numLayers-by-numLayers Boolean matrix
%     [0 1] ... % outputConnect, 1-by-numLayers Boolean vector
%     );
% 
% % number of hidden layer neurons
% neti.layers{1}.size = 1;
% % neti.layers{2}.size = 1;
% 
% % hidden layer transfer function
% neti.layers{1}.transferFcn = 'purelin';
% % neti.layers{2}.transferFcn = 'transig';
% 
% % network training
% neti.trainFcn = 'trainlm';
% neti.performFcn = 'mse';
% 
% neti.trainParam.epochs = 100000;
% 
% neti = configure(neti,inputs,outputs);
% neti = train(neti,inputs,outputs);
% 
% 
% %% READ DATA FROM MOVIES
% 
% clc,clear;tic;
% SLOW = {'Zoom','Zoom0011'; ...
%         'Zoom','Zoom0013'; ...
%         'Zoom','Zoom0014'; ...
%         'Dolly','Dolly0024'; ...
%         'Dolly','Dolly0020'; ...
%         'Dolly','Dolly0009'; ...
%         'Pan','Panning0021'; ...
%         'Pan','Panning0014'; ... 
%         'Pan','Panning0013'; ...
%         'Pedestal','Pedestal0003'; ...
%         'Pedestal','Pedestal0014'; ...
%         'Pedestal','Pedestal0012';  ...
%         'Tilt','Tilt0003'; ...
%         'Tilt','Tilt0016'; ...
%         'Tilt','Tilt0005'; ...
%         'Trucking','Trucking0003'; ...
%         'Trucking','Trucking0011';...
%         'Trucking','Trucking0009'};
% FAST = {'Zoom','Zoom0021'; ...
%         'Zoom','Zoom0019'; ...
%         'Zoom','Zoom0027'; ...
%         'Zoom','Zoom0023'; ...
%         'Dolly','Dolly0011'; ...
%         'Dolly','Dolly0001'; ...
%         'Dolly','Dolly0010'; ...
%         'Dolly','Dolly0012'; ...
%         'Pan','Panning0002'; ...
%         'Pan','Panning0004'; ... 
%         'Pan','Panning0020'; ...
%         'Pedestal','Pedestal0009'; ...
%         'Pedestal','Pedestal0011'; ...
%         'Tilt','Tilt0007'; ...
%         'Tilt','Tilt0011'; ...
%         'Tilt','Tilt0006'; ...
%         'Trucking','Trucking0012'; ...
%         'Trucking','Trucking0013'};
%     
% % Create inputs
% size_input = length(FAST);
% input = [];
% for i = 1:size_input
%     FolderName = strcat('../CAMO_Videos/',FAST{i,1},'/');
%     FileName = strcat(FAST{i,2},'.avi');
%     
%     load(strcat(FolderName,FileName,'_H.mat'));
% 
%     fprintf('%d - %s\n',length(H),FAST{i,2});
%     
%     for i = 1 : 5: length(H)-5
%        
%         currentInput = [];
%         for j = 1:5
%             currentInput = [currentInput ; H(1,1,i+j) ; H(1,2,i+j) ; ...
%                             H(1,3,i+j); H(2,1,i+j) ; H(2,2,i+j) ; ...
%                             H(2,3,i+j); H(3,1,i+j) ; H(3,2,i+j)];
%         end
%         input = [ input , currentInput ]; 
%     end
% end
% 
% size_input = length(SLOW);
% for i = 1:size_input
%     tic;
%     FolderName = strcat('../CAMO_Videos/',SLOW{i,1},'/');
%     FileName = strcat(SLOW{i,2},'.avi');
%     
%     Reader = ReadFunctions;
%     [video,nFrames,vidHeight,vidWidth,~] = Reader.ReadData(FolderName,FileName);
%     mov = Reader.NewMovie(nFrames , vidHeight   ,vidWidth);
%     mov = Reader.ReadMovie(mov , video );
% 
%     fprintf('%d - %s\n',nFrames,SLOW{i,2});
%     
%     H= [];
%     cd('../CameraMotionCode');
%     for n = 1 : nFrames-1
%        
%         im1 = mov(n).cdata;
%         im2 = mov(n+1).cdata;
%         try
%             [T h img] = evalc('homography( im1 , im2 )');
%         catch err
%             h = zeros(3);   
%         end
%         H = cat(3,H,h);
%         toc;
%     end
%     cd('../bin');
%     save(strcat(FolderName,FileName,'_H.mat'),'H');
% end
% 

% %% CREATE INPUT DATA
% 
% FolderName = '../CAMO_Videos/Dolly/';
% FileName = 'Dolly0025.avi';
% load(strcat(FolderName,FileName,'.mat'));
% 
% dollyInput = [];
% for i = 1:5:length(H)-5
%     currentInput = [];
%     for j = 1:5
%         currentInput = [currentInput ; H(1,1,i+j) ; H(1,2,i+j) ; ...
%                         H(1,3,i+j); H(2,1,i+j) ; H(2,2,i+j) ; ...
%                         H(2,3,i+j); H(3,1,i+j) ; H(3,2,i+j)];
%     end
%     dollyInput = [ dollyInput , currentInput ]; 
% end
% 
% dollyMean     =  mean(dollyInput')';
% tiltMean      =  mean(tiltInput')';
% truckingMean  =  mean(truckingInput')';
% zoomMean      =  mean(zoomInput')';
% panMean       =  mean(panInput')';
% pedestalMean  =  mean(pedestalInput')';
% save('../bin/inputs.mat', 'dollyMean' , 'panMean' , ...
%     'pedestalMean' ,'tiltMean' , 'truckingMean' , 'zoomMean' , '-append');
% 
% inputs = [dollyInput(:,1:26) , panInput(:,1:26) , pedestalInput(:,1:26) ...
%           tiltInput(:,1:26) , truckingInput(:,1:26) , zoomInput(:,1:26) ];
% save('../bin/inputs.mat','inputs', 'dollyInput' , 'panInput' , ...
%     'pedestalInput' ,'tiltInput' , 'truckingInput' , 'zoomInput', '-append');

%% TRAIN SVM WITH TWO CAMERA MOTION

clc,clear;
% load('C:\Users\ckocb_000\Thesis\bin\data\H.mat');
load('C:\Users\ckocb_000\Thesis\bin\data\inputs_5.mat');
% input_slow = input_slow(:,randperm(size(input_slow,2)));
% input_fast = input_fast(:,randperm(size(input_fast,2)));

% Create input, output params for SVM
% inputs = [[var(dollyInput(:,1:15))      ;mean(dollyInput(:,1:15))], ...
%           [var(truckingInput(:,1:15))   ;mean(truckingInput(:,1:15))], ...
%           [var(zoomInput(:,1:15))       ;mean(zoomInput(:,1:15))], ...
%           [var(pedestalInput(:,1:15))   ;mean(pedestalInput(:,1:15))],...
%           [var(tiltInput(:,1:15))       ;mean(tiltInput(:,1:15))],...
%           [var(panInput(:,1:15))        ;mean(panInput(:,1:15)) ]]; %first3 = slow = 1, reast = fast = 2
 inputs = [[var(input_slow(:,1:200))  ; mean(input_slow(:,1:200))],...
           [var(input_fast(:,1:200))  ; mean(input_fast(:,1:200))]]';
      
      
% inputs = mapminmax(inputs);
outputs = zeros(length(inputs),1);
for i = 1:length(inputs)
    if      i<200  class = 1;
    else            class = 2;
    end

    outputs(i) = class;
end

% Train
SVMStruct = svmtrain(inputs',outputs,'showplot',true,'kernel_function','rbf');

% Test
input_class1 =  [var(input_slow(:,201:300))  ; mean(input_slow(:,201:300))]';
input_class2 =  [var(input_fast(:,201:300))  ; mean(input_fast(:,201:300))]';    
class1 = svmclassify(SVMStruct,input_class1);
class2 = svmclassify(SVMStruct,input_class2);
size(class1(class1==1))
size(class2(class2==2))


% % Test
% input_class1 = [[var(dollyInput(:,16:23))   ;mean(dollyInput(:,16:23))],...
%                 [var(truckingInput(:,16:23));mean(truckingInput(:,16:23))],...
%                 [var(zoomInput(:,16:23))    ;mean(zoomInput(:,16:23))]]';
% input_class2 = [[var(pedestalInput(:,16:23));mean(pedestalInput(:,16:23))],...
%                 [var(tiltInput(:,16:23))    ;mean(tiltInput(:,16:23))],...
%                 [var(zoomInput(:,16:23))    ;mean(zoomInput(:,16:23))]]';
% % input_class1 = mapminmax(input_class1);
% % input_class2 = mapminmax(input_class2);
% class1 = svmclassify(SVMStruct,input_class1);%, ...
% %                                 dollyInput(9:16,21:29),truckingInput(9:16,21:34),dollyInput(9:16,21:28),...
% %                                 dollyInput(17:24,21:29),truckingInput(17:24,21:34),dollyInput(17:24,21:28),...
% %                                 dollyInput(25:32,21:29),truckingInput(25:32,21:34),dollyInput(25:32,21:28)]');
% class2 = svmclassify(SVMStruct,input_class2);%,...
% %                                 panInput(9:16,21:29),tiltInput(9:16,21:32),zoomInput(9:16,21:30),...
% %                                 panInput(17:24,21:29),tiltInput(17:24,21:32),zoomInput(17:24,21:30),...
% %                                 panInput(25:32,21:29),tiltInput(25:32,21:32),zoomInput(25:32,21:30)]');
% 
% size(class1(class1==1))
% size(class2(class2==2))


