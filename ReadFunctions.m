function F = ReadFunctions
    F.ReadData = @ReadData;
    F.ReadEyeTrackingData = @ReadEyeTrackingData;
    F.ReadMovie = @ReadMovie;
    F.NewMovie = @NewMovie;
    
end


%% PUBLIC FUNCTIONS

function [video,nFrames,vidHeight,vidWidth,vidFPS,vidDuration,frames] = ReadData( foldername,filename )

    % Read video
    video       = VideoReader( strcat(foldername, filename ) );
    nFrames     = video.NumberOfFrames;
    vidHeight   = video.Height;
    vidWidth    = video.Width;
    vidFPS      = video.FrameRate;
    vidDuration = video.Duration;
    frames      = read(video);
    
end

function valueSet = ReadEyeTrackingData(filename)

    % Read eye tracking data
    for subjectNumber = 1:19
        
        if subjectNumber < 10
            currentFileName = strcat('../gaze_hollywood2/samples/00',num2str(subjectNumber),'_' , filename , '.txt');
        else
            currentFileName = strcat('../gaze_hollywood2/samples/0',num2str(subjectNumber),'_' , filename , '.txt');
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
    fid = fopen( '../gaze_hollywood2/geometry.txt' );
    data = textscan( fid , '%f %f %f %d %d' );
    screenResX   = cell2mat( data( 4 ) );
    screenResY   = cell2mat( data( 5 ) );
    fclose( fid );
    
    % Read screen resolution data
    fid = fopen( '../gaze_hollywood2/resolution.txt' );
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

function [mov] = ReadMovie( mov , video )
    nFrames = length(mov);
    for k = 1 : nFrames
        mov(k).cdata = read(video, k);    
    end
end

function [mov] = NewMovie(nFrames , vidHeight   ,vidWidth )
    mov( 1 : floor( nFrames ))= struct('cdata',zeros(vidHeight   ,vidWidth   , 3,'uint8'),'colormap',[]);
end