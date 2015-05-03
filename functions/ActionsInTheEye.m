function F = ActionsInTheEye
    F.ReadEyeTrackingData = @ReadEyeTrackingData;
    F.CalculateMapping = @CalculateMapping;
end

function valueSet = ReadEyeTrackingData(moviePath)

    if regexp( moviePath, '.*ucf.*')
        gazeFolder = 'gaze_ucfsa';
        regexResult = regexp(moviePath,'\','split');
        filename = [regexResult{end-2} '_' regexResult{end-1} '.avi'];
    elseif regexp (moviePath , '.*Hollywood.*' )
        regexResult = regexp(moviePath,'\','split');
        filename = regexResult{end};
        gazeFolder  = 'gaze_hollywood2';
    end

    % Read eye tracking data
    for subjectNumber = 1:19
        
        if subjectNumber < 10
            currentFileName = ['C:\Users\ckocb_000\Downloads\samples/00' num2str(subjectNumber) '_'  filename  '.txt'];
%             currentFileName = strcat('C:\Users\ckocb_000\Downloads/',gazeFolder,'/samples/00',num2str(subjectNumber),'_' , filename , '.txt');
        else
            currentFileName = ['C:\Users\ckocb_000\Downloads\samples/0' num2str(subjectNumber) '_'  filename  '.txt'];
%             currentFileName = strcat('C:\Users\ckocb_000\Downloads/',gazeFolder,'/samples/0',num2str(subjectNumber),'_' , filename , '.txt');
        end
        
        if  ~exist(currentFileName, 'file')
            continue;
        end
        
        fid = fopen( currentFileName );
        data = textscan( fid , '%d %f %f %f %f %s' );
        fclose( fid );
        
        event = cell2mat(data{:,6});
        event = event~='B';
        
        currentTimeStamp = cell2mat( data( 1 ) );
        currentXScreen = cell2mat( data( 4 ) );
        currentYScreen = cell2mat( data( 5 ) );
        timeStamp{subjectNumber} = currentTimeStamp(event);
        xScreen{subjectNumber} = currentXScreen(event);
        yScreen{subjectNumber} = currentYScreen(event);
    
    end
    
    % Read calibration of the experiment setup
%     fid = fopen( strcat('F:/Thesis/',gazeFolder,'/geometry.txt') );
    fid = fopen( strcat('C:\Users\ckocb_000\Downloads/geometry.txt') );
    data = textscan( fid , '%f %f %f %d %d' );
    screenResX   = cell2mat( data( 4 ) );
    screenResY   = cell2mat( data( 5 ) );
    fclose( fid );
    
    % Read screen resolution data
%     fid = fopen( strcat('F:/Thesis/',gazeFolder,'/resolution.txt' ));
    fid = fopen( strcat('C:\Users\ckocb_000\Downloads/resolution.txt') );
    data = textscan( fid , '%s %d %d %f' );
    fclose( fid );
    vidNames =  [ data{ 1 } ] ;
    index = (strcmp( filename , vidNames ) );
    videoResX   = cell2mat( data( 2 ) );
    videoResX   = videoResX( find( index == 1 ) );
    videoResY   = cell2mat( data( 3 ) );
    videoResY   = videoResY( find( index == 1 ) );
   
    valueSet = struct();     
    valueSet.timeStamp = timeStamp;
    valueSet.xScreen = xScreen;
    valueSet.yScreen = yScreen;
    valueSet.screenResX = screenResX;
    valueSet.screenResY = screenResY;
    valueSet.videoResX = videoResX;
    valueSet.videoResY = videoResY;
    
end

function eyes = CalculateMapping(mappingStruct)

    videoResX  = mappingStruct.videoResX;
    videoResY  = mappingStruct.videoResY;
    screenResX = mappingStruct.screenResX;
    screenResY = mappingStruct.screenResY;
    xScreen    = mappingStruct.xScreen;
    yScreen    = mappingStruct.yScreen;
    timeStamp  = mappingStruct.timeStamp;
    vidDuration= mappingStruct.vidDuration;
    nFrames    = mappingStruct.nFrames;

    a           = double( videoResX ) / double( screenResX ) ; 
    b           = ( double( screenResY ) - double( videoResY ) / a ) / 2;
    eyes = [];
    for subjectNumber = 1:19
        xHeight     = a * xScreen{subjectNumber};
        yHeight     = a * ( yScreen{subjectNumber} - b );
        frames      = int32(round( nFrames * double(timeStamp{subjectNumber}) / ( vidDuration * 1000000 ) )) +1;
        eyes        = [ eyes; repmat(subjectNumber,size(frames)) frames yHeight xHeight];
    end

    eyes = eyes( eyes(:,2) <= nFrames & ...
                 eyes(:,3) > 1 & ...
                 eyes(:,3) < videoResY & ...
                 eyes(:,4) > 1 & ...
                 eyes(:,4) < videoResX,: );
%     eyes = sortrows(eyes,1);

end