function F = SaliencyFunctions
%     TimerScript;
    F.SaveEyeFixationForLeMeur = @SaveEyeFixationForLeMeur;
    F.SaveFramesForLeMeur = @SaveFramesForLeMeur;
    F.CalculateSaliency = @CalculateSaliency;
    F.PlotROC = @PlotROC;
    F.ReadRValuesFromFiles = @ReadRValuesFromFiles;
    F.SaveOpticalFlowData = @SaveOpticalFlowData;
    F.TrainGaussCenterPredicter = @TrainGaussCenterPredicter;
    F.SaveHomographyData = @SaveHomographyData;
end

%% PUBLIC FUNCTIONS
function SaveEyeFixationForLeMeur(folderName, fileName)

    movieNames = {'actioncliptrain00240.avi', ...
                'actioncliptrain00264.avi', ...
                'actioncliptrain00135.avi', ...
                'actioncliptrain00141.avi', ...
                'actioncliptrain00172.avi', ...
                'actioncliptrain00234.avi', ...
                'actioncliptrain00121.avi', ...
                'actioncliptrain00094.avi', ...
                'actioncliptrain00073.avi', ...
                'actioncliptrain00039.avi', ...
                'actioncliptrain00057.avi', ...
                'actioncliptrain00681.avi', ...
                'actioncliptrain00218.avi', ...
                'actioncliptrain00088.avi', ...
                'actioncliptrain00134.avi'};

    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    
    Reader = ReadFunctions;
    Util = UtilFunctions;
    a = 1;

    for i = 1:size(movieNames,2)
        
        fileName = movieNames{i};
        
        [~,nFrames,vidHeight,vidWidth,~,vidDuration] = Reader.ReadData(folderName,fileName);
        resultMap = Reader.ReadEyeTrackingData(fileName);
        resultMap.vidDuration = vidDuration;
        resultMap.nFrames = nFrames;
        saliencyPoints = Util.CalculateMapping(resultMap);

        for k = 1:nFrames
            
            strToWrite = [];
            for s = 1:19
                indices = find(saliencyPoints( : , 2 ) == k & saliencyPoints( : , 1 ) == s);

                for t = 1:size(indices,1)
                    x = saliencyPoints(indices(t),4);
                    y = saliencyPoints(indices(t),3);
                    strToWrite = [strToWrite num2str(x) ' ' num2str(y) ' 20 '] ;
                end
                if ~isempty(indices)
                    strToWrite = [strToWrite '-1 -1 -1\n'];
                end

            end

            strToWrite = [strToWrite '-1 -1 -1'];

            mkdir('..\forLeMeur\massVideos')
            fid = fopen( strcat('..\forLeMeur\massVideos\frame_',num2str(a),'.stat'), 'wt' );
            fprintf(fid,strToWrite);
            fclose(fid);

            a=a+1;
              
        end
        
    end
    
end

function SaveFramesForLeMeur(folderName, fileName)

    movieNames = {'actioncliptrain00240.avi', ...
                'actioncliptrain00264.avi', ...
                'actioncliptrain00135.avi', ...
                'actioncliptrain00141.avi', ...
                'actioncliptrain00172.avi', ...
                'actioncliptrain00234.avi', ...
                'actioncliptrain00121.avi', ...
                'actioncliptrain00094.avi', ...
                'actioncliptrain00073.avi', ...
                'actioncliptrain00039.avi', ...
                'actioncliptrain00057.avi', ...
                'actioncliptrain00681.avi', ...
                'actioncliptrain00218.avi', ...
                'actioncliptrain00088.avi', ...
                'actioncliptrain00134.avi'};    

    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    
    a = 1;
    Reader = ReadFunctions;
    
    for i = 1:size(movieNames,2)
        
        fileName = movieNames{i};
       
        [~,nFrames,~,~,~,~,frames] = Reader.ReadData(folderName,fileName);
        
    %     mkdir(strcat('..\forLeMeur\',fileName));
        mkdir('..\forLeMeur\massVideos');
        for k = 1:nFrames
    %         imwrite( frames(:,:,:,k), ...
    %             strcat('..\forLeMeur\',fileName,'\frame_',num2str(k),'.bmp') , 'bmp');
            imwrite( frames(:,:,:,k), ...
                        strcat('..\forLeMeur\massVideos\frame_',num2str(a),'.bmp') , 'bmp');
            a=a+1;
        end
        
    end
    
end

function CalculateSaliency(folderName, mode)

    if nargin<2
        mode = 'gbvs';
    end
    
    movieNames = {'actioncliptrain00240.avi', ...
                'actioncliptrain00264.avi', ...
                'actioncliptrain00135.avi', ...
                'actioncliptrain00141.avi', ...
                'actioncliptrain00172.avi', ...
                'actioncliptrain00234.avi', ...
                'actioncliptrain00121.avi', ...
                'actioncliptrain00094.avi', ...
                'actioncliptrain00073.avi', ...
                'actioncliptrain00039.avi', ...
                'actioncliptrain00057.avi', ...
                'actioncliptrain00681.avi', ...
                'actioncliptrain00218.avi', ...
                'actioncliptrain00088.avi', ...
                'actioncliptrain00134.avi'};    

    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    
    Reader = ReadFunctions;
    Video = VideoFunctions;

    mkdir(strcat('..\forLeMeur\massVideos\CalculatedSaliency_',mode));
       
    tic;   
    a = 1;
    for i = 1:size(movieNames,2)
        
        fileName = movieNames{i};

        [video,nFrames,vidHeight,vidWidth,~,~,frames] = Reader.ReadData(folderName,fileName);
        mov = Reader.NewMovie(nFrames , vidHeight   ,vidWidth);
        mov = Reader.ReadMovie(mov , video );

        if strcmp( mode , 'gauss' )
            load(strcat('data/OpticalFlowMaps',num2str(i),'.mat'));
            
            for k = 1:nFrames

                [gaussCenterX,gaussCenterY] = GetGaussCenter(opt(:,:,k,1),opt(:,:,k,2));
                gaussCenterX = gaussCenterX * vidWidth;
                gaussCenterY = gaussCenterY * vidHeight;
                
                offsetX = gaussCenterX - vidWidth/2;
                offsetY = gaussCenterY - vidHeight/2;
                
                currentGauss = MyGauss(vidHeight,vidWidth,vidHeight/5,offsetX,offsetY);
                videoSaliencyMap = mat2gray(currentGauss);

                imwrite( videoSaliencyMap, ...
                    strcat('..\forLeMeur\massVideos\CalculatedSaliency_',mode,'\frame_',num2str(a),'.bmp') , 'bmp');
                a = a + 1;
            end
            
        elseif strcmp( mode , 'cigdem' ) 
            
%             avgSaliency = Video.CalculateMeanSaliency( nFrames , saliencyPoints ); 
%             shotBoundaries = Video.DetectShotBoundaries( mov );
%             load(strcat('data/OpticalFlowMaps',num2str(i),'.mat'));
%             [avgFlowOptical] = Video.CreateFlow( 'optical' , shotBoundaries , avgSaliency , mov , 5, opt ); 

            load(strcat('data/OpticalFlowMaps',num2str(i),'.mat'));
            [~ , saliencyMap] = ReadStaticSaliency(a,nFrames,'judd');
            
            for k = 1:nFrames
                
                [gaussCenterX,gaussCenterY] = GetGaussCenter(opt(:,:,k,1),opt(:,:,k,2));
                gaussCenterX = gaussCenterX * vidWidth;
                gaussCenterY = gaussCenterY * vidHeight;
                
                offsetX = gaussCenterX - vidWidth/2;
                offsetY = gaussCenterY - vidHeight/2;
                
                currentGauss = MyGauss(vidHeight,vidWidth,vidHeight/3,offsetX,offsetY);

                videoSaliencyMap = mat2gray(saliencyMap(:,:,k).*currentGauss);

                imwrite( videoSaliencyMap, ...
                    strcat('..\forLeMeur\massVideos\CalculatedSaliency_',mode,'\frame_',num2str(a),'.bmp') , 'bmp');
                a = a + 1;
            end

        else

            [~ , saliencyMap] = Video.CalculateStaticSaliency(mov,mode);

            for k = 1:nFrames
                
                videoSaliencyMap = saliencyMap(:,:,k);
                imwrite( videoSaliencyMap, ...
                    strcat('..\forLeMeur\massVideos\CalculatedSaliency_',mode,'\frame_',num2str(a),'.bmp') , 'bmp');
                a = a + 1;
                
            end

        end
        toc;
    end
    
end

function [score,R] = PlotROC(algorithm,plotColor)

    if nargin < 2
        plotColor = 'b';
    end
    if nargin == 0
        algorithm = 'sr';
    end
    
    Evaluate = EvalFunctions; 
    a = 1;
    for k = [397:597,607:824,862:924,947:1041,1072:1185,1213:1346,1502:1677,1678:1862,1871:1906,1937:2211,2246:2465,2469:2751]
        
        eyeFolder = strcat('../forLeMeur/massVideos/leMeurFixationMaps/frame_',num2str(k),'_saliencyMapGT.bmp');
        salFolder = strcat('../forLeMeur/massVideos/CalculatedSaliency_',algorithm,'/frame_',num2str(k),'.bmp');
        
        currentEye = imread(eyeFolder);
        currentSaliency = imread(salFolder);

        eyeFixation{a} = currentEye(:,:,1);
        saliency{a} = currentSaliency(:,:,1);
        a = a + 1;

    end
    
    [score,R] = Evaluate.CalculateAUCscoreVideo( saliency, eyeFixation, algorithm );
    
    plot(R(:,1),R(:,2),plotColor);
    
end

function outputImage = MyHistMatch(inputImage, targetImage)

    [c,x] = imhist(targetImage); 
    outputImage = histoMatch(inputImage, c, x);
     
end

function f = MyGauss(vidHeight, vidWidth, sigma, centerx , centery)

    [x y]=meshgrid( round(-vidWidth/2)+1:round(vidWidth/2), ...
                    round(-vidHeight/2)+1:round(vidHeight/2));
                
    f = exp( ( -(x-centerx).^2 / (2*sigma^2) ) + ...
           ( -(y-centery).^2 / (2*sigma^2) ) )  ;
    f = f./sum(f(:));
    
end

function [saliencyPoints , saliencyMap] = ReadStaticSaliency(a,nFrames,mode)

    for k = 1:nFrames
        saliencyMap(:,:,k) = im2double(imread( strcat('..\forLeMeur\massVideos\CalculatedSaliency_',mode,'\frame_',num2str(a),'.bmp')));
        a = a + 1;
    end
    
    saliencyPoints = [];
    for k = 1 : nFrames
        frame_map = saliencyMap(:,:,k);
        [r c] = find(frame_map>0);
        sizer = size(r);
        saliencyPoints = [saliencyPoints ; [repmat(k,sizer,1) r c]];

    end
    

end

function [score,totalR] = ReadRValuesFromFiles(algorithm)

    totalR = zeros(12,2);
    for k = 100:100:2000
        R = [];
        load(strcat('data\R_',algorithm,'_',num2str(k),'.mat'));
        totalR = totalR + R;
    end
    
    totalR = totalR / 2000;
    score = trapz(flipdim(totalR(:,1),1),flipdim(totalR(:,2),1));
    
end

function SaveOpticalFlowData
    movieNames = {'actioncliptrain00240.avi', ...
                'actioncliptrain00264.avi', ...
                'actioncliptrain00135.avi', ...
                'actioncliptrain00141.avi', ...
                'actioncliptrain00172.avi', ...
                'actioncliptrain00234.avi', ...
                'actioncliptrain00121.avi', ...
                'actioncliptrain00094.avi', ...
                'actioncliptrain00073.avi', ...
                'actioncliptrain00039.avi', ...
                'actioncliptrain00057.avi', ...
                'actioncliptrain00681.avi', ...
                'actioncliptrain00218.avi', ...
                'actioncliptrain00088.avi', ...
                'actioncliptrain00134.avi'};   
    
    Reader = ReadFunctions;
    Video = VideoFunctions;
   
    for i = 1:size(movieNames,2)
        
        fileName = movieNames{i};

        [video,nFrames,vidHeight,vidWidth,~,~,frames] = Reader.ReadData('../Hollywood2-actions/Hollywood2/AVIClips/',fileName);
        mov = Reader.NewMovie(nFrames , vidHeight   ,vidWidth);
        mov = Reader.ReadMovie(mov , video );
        
        [flowX,flowY] = Video.OpticalFlowMap( mov ); 
        opt = cat(4,flowX,flowY);
        save(strcat('data/OpticalFlowMaps',num2str(i),'.mat'),'opt');
        
    end
    
end

function SaveHomographyData
    movieNames = {'actioncliptrain00240.avi', ...
                'actioncliptrain00264.avi', ...
                'actioncliptrain00135.avi', ...
                'actioncliptrain00141.avi', ...
                'actioncliptrain00172.avi', ...
                'actioncliptrain00234.avi', ...
                'actioncliptrain00121.avi', ...
                'actioncliptrain00094.avi', ...
                'actioncliptrain00073.avi', ...
                'actioncliptrain00039.avi', ...
                'actioncliptrain00057.avi', ...
                'actioncliptrain00681.avi', ...
                'actioncliptrain00218.avi', ...
                'actioncliptrain00088.avi', ...
                'actioncliptrain00134.avi'};   
    
    Reader = ReadFunctions;
    Video = VideoFunctions;
   
    for i = 1:size(movieNames,2)
        
        fileName = movieNames{i};

        [video,nFrames,vidHeight,vidWidth,~,~,frames] = Reader.ReadData('../Hollywood2-actions/Hollywood2/AVIClips/',fileName);
        mov = Reader.NewMovie(nFrames , vidHeight   ,vidWidth);
        mov = Reader.ReadMovie(mov , video );
        
        h = Video.CalculateHomography( mov ); 
        save(strcat('data/HomographyMatrices',num2str(i),'.mat'),'h');
        
    end
    
end

function TrainGaussCenterPredicter()

    horPatchNumber = 5;
    verPatchNumber = 3;

    inputsX = [];    outputsX = [];    inputsTestX = [];    outputsTestX = [];
    inputsY = [];    outputsY = [];    inputsTestY = [];    outputsTestY = [];
    
    saliencyMapIndex = 1;
    for k = 1:15
        
        load(strcat('data/OpticalFlowMaps',num2str(k),'.mat'));
        [frameHeight , frameWidth, ~ , ~] = size(opt);
        verPatchSize = floor(frameHeight / verPatchNumber);
        horPatchSize = floor(frameWidth / horPatchNumber);
        
        for t = 2:size(opt,3)-1
            
            eyeFolder = strcat('../forLeMeur/massVideos/leMeurFixationMaps/frame_',num2str(saliencyMapIndex),'_saliencyMapGT.bmp');
            
            currentEye = imread(eyeFolder);
            currentEye = currentEye(:,:,1);
            [~,ind] = max(currentEye(:));
            [y,x]=ind2sub(size(currentEye),ind);
%             y = ceil(y/verPatchSize);
%             x = ceil(x/horPatchSize);
            
            opticalFlowX = opt(:,:,t-1:t+1,1);
            opticalFlowY = opt(:,:,t-1:t+1,2);
            opticalFlowX = mat2gray(opticalFlowX(:,:,:));
            opticalFlowY = mat2gray(opticalFlowY(:,:,:));
%             currentOpticalFlow = sqrt(opticalFlowX.^2 + opticalFlowY.^2);
            
            patchIndex = 1;
            for i = 1:horPatchSize:frameWidth-10
                for j = 1:verPatchSize:frameHeight-10
                    currentPatchX = opticalFlowX(j:j+verPatchSize-1,i:i+horPatchSize-1,:);
                    currentPatchY = opticalFlowY(j:j+verPatchSize-1,i:i+horPatchSize-1,:);
                    patchMeanX(patchIndex) = (mean(currentPatchX(:)));
                    patchMeanY(patchIndex) = (mean(currentPatchY(:)));
                    patchIndex = patchIndex + 1;
                end
            end

            if sum(patchMeanX) ~= 0
                
                % normalize the input/output data
%                 patchMeanX =  patchMeanX/max(patchMeanX);
%                 patchMeanY =  patchMeanY/max(patchMeanY);
                patchMeanX =  mapminmax(patchMeanX);
                patchMeanY =  mapminmax(patchMeanY);
                x =  x/frameWidth;
                y =  y/frameHeight;

                if saliencyMapIndex ~= [397:597,607:824,862:924,947:1041,1072:1185,1213:1346,1502:1677,1678:1862,1871:1906,1937:2211,2246:2465,2469:2751]
                    inputsX = [inputsX,patchMeanX'];
                    inputsY = [inputsY,patchMeanY'];
                    outputsX = [outputsX, x];
                    outputsY = [outputsY, y];
                else
                    inputsTestX = [inputsTestX,patchMeanX'];
                    inputsTestY = [inputsTestY,patchMeanY'];
                    outputsTestX = [outputsTestX, x];
                    outputsTestY = [outputsTestY, y];
                end
                
            end
            
            saliencyMapIndex=saliencyMapIndex+1;

        end

    end
    
%     save('inputs.mat','inputsX','inputsTestX','inputsY','inputsTestY', ...
%                       'outputsX','outputsTestX','outputsY','outputsTestY');
%     load('inputs_minusplus.mat');
   
    nnetX = feedforwardnet(20);
%     nnet.trainFcn = 'trainscg';
    [nnetX,tr] = train(nnetX,inputsX,outputsX); 
    predictedClassX = nnetX(inputsTestX);    
    predictedClassX = predictedClassX';
    actualClass = outputsTestX';
    errorRateX = mean(abs(predictedClassX - actualClass))
    
    nnetY = feedforwardnet(20);
%     nnet.trainFcn = 'trainscg';
    [nnetY,tr] = train(nnetY,inputsY,outputsY); 
    predictedClassY = nnetY(inputsTestY);    
    predictedClassY = predictedClassY';
    actualClass = outputsTestY';
    errorRateY = mean(abs(predictedClassY - actualClass))

%     save('data/GaussCenterPredicter.mat' , 'nnetX' , 'nnetY');
    
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









