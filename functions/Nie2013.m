function F = Nie2013
    F.WriteMaskFile = @WriteMaskFile;
    F.GridToOff = @GridToOff;
end

function WriteMaskFile(data,outputPath)

    fid = fopen([outputPath '\mask.txt'],'wt');
%     fid = fopen('C:\Users\ckocb_000\Documents\Visual Studio 2012\Projects\TetGen\Debug\mask.txt','wt');
    fprintf('Group Manager\n');
    fprintf(fid, '%d\n', data.NumberOfGroups);
    
    for i = 0:data.NumberOfGroups-1
        
        currentGroup = data.Groups{i+1};
        
        if isempty(currentGroup)
            fprintf(fid, 'Group\n%d\nNone\n', i-1);
            continue;
        end
        
        currentFrameNum = currentGroup.NumberOfFrames;
        fprintf(fid, 'Group\n%d\nGraphicsGroup\n%d', i-1, currentFrameNum);  
        
        for j = 0:currentFrameNum-1
            
            currentFrame = currentGroup.Frames{j+1};
            
            if isempty(currentFrame)
                fprintf(fid, 'Frame\n%d\nNone\n', j-1);
                continue;
            end
            
            for k = 1:size(currentFrame,1)
                fprintf(fid, '%d %d\n', currentFrame(k,1),currentFrame(k,2));  
            end
            
        end
        
    end
    fclose(fid);
end

function GridToOff
    % Create a grid
    [x,y,z] = meshgrid(1:10:100,1:10:300,1:10:40);
    vertex = [x(:) y(:) z(:)]';
    face = delaunayTriangulation(x(:), y(:), z(:)); 
    
    fid = fopen('C:\Users\ckocb_000\Documents\Visual Studio 2012\Projects\TetGen\Debug\test.node','wt');
    fprintf(fid, '%d 3 0 0\n', size(vertex,2));
    for i = 1:size(vertex,2)
        fprintf(fid, '%d %f %f %f\n', i, vertex(1,i),vertex(2,i),vertex(3,i));
    end
    fclose(fid);
end