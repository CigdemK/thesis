function F = ReadFunctions
    F.ReadData = @ReadData;
    F.ReadEyeTrackingData = @ReadEyeTrackingData;
    F.ReadMovie = @ReadMovie;
    F.NewMovie = @NewMovie;
    
end


%% PUBLIC FUNCTIONS

function [video,nFrames,vidHeight,vidWidth,vidFPS] = ReadData( foldername,filename )

    % Read video
    video       = VideoReader( strcat(foldername, filename ) );
    nFrames     = video.NumberOfFrames;
    vidHeight   = video.Height;
    vidWidth    = video.Width;
    vidFPS      = video.FrameRate;
    
end

function ReadEyeTrackingData(foldername,filename)

    load('Data.mat');
    % Read eye tracking data
    fid = fopen( strcat(foldername,'006_' , filename , '.txt') );
    data = textscan( fid , '%d %f %f %f %f %s' );
    fclose( fid );
    timeStamp   = cell2mat( data( 1 ) );
    xScreen     = cell2mat( data( 4 ) );
    yScreen     = cell2mat( data( 5 ) );
    
    % Read calibration of the experiment setup
    fid = fopen( '../holywood2/geometry.txt' );
    data = textscan( fid , '%f %f %f %d %d' );
    screenResX   = cell2mat( data( 4 ) );
    screenResY   = cell2mat( data( 5 ) );
    fclose( fid );
    
    % Read screen resolution data
    fid = fopen( '../holywood2/resolution.txt' );
    data = textscan( fid , '%s %d %d %f' );
    fclose( fid );
    vidNames =  [ data{ 1 } ] ;
    index = (strcmp( filename , vidNames ) );
    videoResX   = cell2mat( data( 2 ) );
    videoResX   = videoResX( find( index == 1 ) );
    videoResY   = cell2mat( data( 3 ) );
    videoResY   = videoResY( find( index == 1 ) );
   
    save('Data.mat');  
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