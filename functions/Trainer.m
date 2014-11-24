function F = Trainer
    F.TrainGaussCenterPredicter = @TrainGaussCenterPredicter;
    F.PrepareIOForNN_OpticalFlow = @PrepareIOForNN_OpticalFlow;
    F.PrepareIOForNN_Trajectory = @PrepareIOForNN_Trajectory;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Public Functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [inputsX,inputsY,outputsX,outputsY] = PrepareIOForNN_OpticalFlow(folderName,movieNames)

    inputsX = [];    outputsX = [];    
    inputsY = [];    outputsY = [];   
    
    horPatchNumber = 5;
    verPatchNumber = 3;
    
    tic;
    for k = 1:size(movieNames,2)
        
        load(strcat(folderName,movieNames{k},'_opticalFlow.mat'));
        [frameHeight , frameWidth, frameNumber , ~] = size(opt);
        verPatchSize = floor(frameHeight / verPatchNumber);
        horPatchSize = floor(frameWidth / horPatchNumber);
        
        for t = 2:frameNumber-1
            
            regexResult = regexp(movieNames(k),'/','split');
            eyeFolder = strcat(folderName,regexResult{1}{1},'/', ...
                regexResult{1}{2},'/frame_',num2str(t),'_saliencyMapGT.bmp');
            
            currentEye = imread(eyeFolder);
            currentEye = currentEye(:,:,1);
            [~,ind] = max(currentEye(:));
            [y,x]=ind2sub(size(currentEye),ind);
%             y = ceil(y/verPatchSize);
%             x = ceil(x/horPatchSize);
            
            opticalFlowX = opt(:,:,t,1);
            opticalFlowY = opt(:,:,t,2);
%             opticalFlowX = mat2gray(opticalFlowX(:,:,:));
%             opticalFlowY = mat2gray(opticalFlowY(:,:,:));
            
            patchIndex = 1;
            for i = 1:horPatchSize:frameWidth-10
                for j = 1:verPatchSize:frameHeight-10
                    currentPatchX = opticalFlowX(j:j+verPatchSize-1,i:i+horPatchSize-1);
                    currentPatchY = opticalFlowY(j:j+verPatchSize-1,i:i+horPatchSize-1);
                    patchMeanX(patchIndex) = (mean(currentPatchX(:)));
                    patchMeanY(patchIndex) = (mean(currentPatchY(:)));
                    patchIndex = patchIndex + 1;
                end
            end


            % normalize the input/output data
            patchMeanX =  patchMeanX/max(patchMeanX);
            patchMeanY =  patchMeanY/max(patchMeanY);
%             patchMeanX =  mapminmax(patchMeanX);
%             patchMeanY =  mapminmax(patchMeanY);
%             max(patchMeanX)
%             min(patchMeanX)
            x =  x/frameWidth;
            y =  y/frameHeight;

            inputsX = [inputsX,patchMeanX'];
            inputsY = [inputsY,patchMeanY'];
            outputsX = [outputsX, x];
            outputsY = [outputsY, y];

        end
        toc;
    end
end

function [inputs,outputs] = PrepareIOForNN_Trajectory(folderName,movieNames)

    inputs = [];    
    outputs = [];    
    
    for k = 1:size(movieNames,2)
        
        load(strcat(folderName,movieNames{k},'_ImprovedTrajectory_TrajectoriesByFrame.mat'));
        load(strcat(folderName,movieNames{k},'_opticalFlow.mat'));
        [frameHeight , frameWidth, frameNumber , ~] = size(opt);
         
        for t = 1:frameNumber-1
            
            regexResult = regexp(movieNames(k),'/','split');
            eyeFolder = strcat(folderName,regexResult{1}{1},'/', ...
                regexResult{1}{2},'/frame_',num2str(t),'_saliencyMapGT.bmp');
            
            currentEye = imread(eyeFolder);
            currentEye = currentEye(:,:,1);
            [~,ind] = max(currentEye(:));
            [y,x]=ind2sub(size(currentEye),ind);
            
            meanTrajectoriesX   = mean(trajectoriesByFrames{t}(1,:));
            meanTrajectoriesY   = mean(trajectoriesByFrames{t}(2,:));
%             meanTrajectoriesHOF = mean(trajectoriesByFrames{t}(3:end,:),2);
            
            % normalize the input/output data
            x =  x/frameWidth;
            y =  y/frameHeight;
            meanTrajectoriesX =  meanTrajectoriesX/frameWidth;
            meanTrajectoriesY =  meanTrajectoriesY/frameHeight;

%             inputs = [inputs,[meanTrajectoriesX;meanTrajectoriesY;meanTrajectoriesHOF]];
            inputs = [inputs,[meanTrajectoriesX;meanTrajectoriesY]];            
            outputs = [outputs, [x;y]];

        end
    end
    save('inputs.mat','inputs','outputs');
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Private Functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TrainGaussCenterPredicter()

%     [inputsX,inputsY,outputsX,outputsY] = PrepareIOForNN();
    load('inputs.mat')

    nnet = feedforwardnet(10);
%     nnet.trainFcn = 'trainscg';
    [nnet,tr] = train(nnet,inputs(:,1:400),outputs(:,1:400)); 
    save('data/GaussCenterPredicter.mat' , 'nnet');
    
%     predictedClass = nnet(inputs(:,1:400));    
%     predictedClass = predictedClass';
%     actualClass = outputs(:,1:400)';
%     errorRate = mean(abs(predictedClass - actualClass))
%     
%     predictedClass = nnet(inputs(:,400:end));    
%     predictedClass = predictedClass';
%     actualClass = outputs(:,400:end)';
%     errorRate = mean(abs(predictedClass - actualClass))
    
end

function [x,y] = GetGaussCenter(opticalFlowX,opticalFlowY)

    [frameHeight , frameWidth] = size(opticalFlowX);
    verPatchSize = floor(frameHeight / 3);
    horPatchSize = floor(frameWidth / 5);
    
%     currentOpticalFlow = sqrt(opticalFlowX.^2 + opticalFlowY.^2);

    patchIndex = 1;
    for i = 1:horPatchSize:frameWidth-10
        for j = 1:verPatchSize:frameHeight-10
            currentPatchX = opticalFlowX(j:j+verPatchSize-1,i:i+horPatchSize-1);
            currentPatchY = opticalFlowY(j:j+verPatchSize-1,i:i+horPatchSize-1);
            patchesX(patchIndex) = mean(currentPatchX(:));
            patchesY(patchIndex) = mean(currentPatchY(:));
            patchIndex = patchIndex+1;
        end
    end

    load('data/GaussCenterPredicter.mat');
    x = nnetX(patchesX');
    y = nnetY(patchesY');

end