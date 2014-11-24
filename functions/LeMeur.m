function F = LeMeur 
    F.SaveEyeFixation = @SaveEyeFixation;
    F.SaveFrames = @SaveFrames;
end

function SaveEyeFixation(moviePath)

    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    
    video = VideoReader( moviePath );
    
    AITE = ActionsInTheEye;
    resultMap = AITE.ReadEyeTrackingData(moviePath);
    resultMap.vidDuration = video.Duration;
    resultMap.nFrames = video.NumberOfFrames;
    saliencyPoints = AITE.CalculateMapping(resultMap);

    for k = 1:video.NumberOfFrames

        strToWrite = [];
        for subject = 1:19
            indices = find(saliencyPoints( : , 2 ) == k & saliencyPoints( : , 1 ) == subject);

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

        regexResult = regexp(moviePath,'\','split');
        fid = fopen( [regexResult{1} '\' regexResult{2} '\' regexResult{3} '\' ...
            regexResult{4} '\' regexResult{5} '\' regexResult{6} '\frame_',num2str(k),'.stat'], 'wt' );
        fprintf(fid,strToWrite);
        fclose(fid);

    end
end

function SaveFrames(moviePath)

    warning('off', 'MATLAB:MKDIR:DirectoryExists');

    video = VideoReader( moviePath );
    frames = read(video);

    regexResult = regexp(moviePath,'\','split');
    folderName = [regexResult{1} '\' regexResult{2} '\' regexResult{3} '\' ...
           regexResult{4} '\' regexResult{5} '\' regexResult{6}];
    for k = 1:video.NumberOfFrames
        imwrite( frames(:,:,:,k), ...
            [folderName '\frame_' num2str(k) '.png'] , 'png');
    end

end
