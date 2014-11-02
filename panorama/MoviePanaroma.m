function [panMov] = MoviePanaroma(FolderName , FileName )
    
    % Read data
    Reader = ReadFunctions;
    [video,nFrames,vidHeight,vidWidth,~] = Reader.ReadData(FolderName,FileName);
    mov = Reader.NewMovie(nFrames , vidHeight   ,vidWidth);
    mov = Reader.ReadMovie(mov , video );
    
    Cropper = CropFunctions;
    mov = Cropper.RemoveBlackBars(mov);
    
    % Create panaroma
    tic;panMov = Panorama(mov(1).cdata,mov(2).cdata);toc;
    for n = 3:20
        panMov = Panorama(panMov,mov(n).cdata); toc;
    end
    
end