function F = ImprovedTrajectories
    F.RunImpTrajectories = @RunImpTrajectories;
    F.AnalyzeOutput = @AnalyzeOutput;
    F.GetTrajectoriesByFrame = @GetTrajectoriesByFrame;
end

function [status,output] = RunImpTrajectories(videoPath)
    
    exePath = 'C:\Users\ckocb_000\Documents\Visual Studio 2010\Projects\ImprovedTrajectories\Debug\ImprovedTrajectories.exe';
    [status,output] = system(['"',exePath,'" "',videoPath,'"']);
    
end

function trajectories = AnalyzeOutput(exeOutput)
      
    tempTrajectories = regexp(exeOutput,'\n','split');
    numberOfTrajectories = size(tempTrajectories,2);
    
    U = UtilFunctions;
    trajectories = cell(numberOfTrajectories-2,1);
    for k = 1:numberOfTrajectories-2 % the last element is empty, and the previous one contains info about the movie

        currentTempTrajectory = regexp(tempTrajectories{k},'\s','split');
        currentTempTrajectory = U.flattenCellArray(currentTempTrajectory);
        currentTempTrajectory = str2double(currentTempTrajectory);     
        
        currentTrajectory = struct();
        
        currentTrajectory.frameNum  = currentTempTrajectory(1);
        currentTrajectory.mean_x    = currentTempTrajectory(2);   
        currentTrajectory.mean_y    = currentTempTrajectory(3);    
        currentTrajectory.var_x     = currentTempTrajectory(4); 
        currentTrajectory.var_y     = currentTempTrajectory(5); 
        currentTrajectory.length    = currentTempTrajectory(6); 
        currentTrajectory.scale     = currentTempTrajectory(7);     
        currentTrajectory.x_pos     = currentTempTrajectory(8);       
        currentTrajectory.y_pos     = currentTempTrajectory(9);       
        currentTrajectory.t_pos     = currentTempTrajectory(10);  
        
        trajectoryLenght =size( currentTempTrajectory, 2 ) - 10 - 96 - 108 - 96 - 96; %length in the array = 2*actualLength
        if trajectoryLenght > 1
            currentTrajectory.trajectory= reshape(  currentTempTrajectory( 11 : 10+trajectoryLenght ) , [2 trajectoryLenght/2] ); 
        else
            currentTrajectory.trajectory = [];
        end
%         currentTrajectory.HOG       = reshape(  currentTempTrajectory( trajectoryLenght+11  : trajectoryLenght+106 ) , [8 2 2 3] );         
        currentTrajectory.HOF       = reshape(  currentTempTrajectory( trajectoryLenght+107 : trajectoryLenght+214 ) , [9 2 2 3] );          
%         currentTrajectory.MBHx      = reshape(  currentTempTrajectory( trajectoryLenght+215 : trajectoryLenght+310 ) , [8 2 2 3] );
%         currentTrajectory.MBHy      = reshape(  currentTempTrajectory( trajectoryLenght+311 : trajectoryLenght+406 ) , [8 2 2 3] );
        
        trajectories{k} = currentTrajectory;
        
    end
    
end

function trajectoriesByFrames = GetTrajectoriesByFrame(trajectories)

    numberOfTrajectories = size(trajectories,1);
    numberOfFrames = trajectories{end}.frameNum;
    trajectoriesByFrames = cell(numberOfFrames,1);
    
    for i = 1: numberOfTrajectories

        currentTrajectory = trajectories{i};
        trajectoryLength = size(currentTrajectory.trajectory,2);
        trajectoryEnd = currentTrajectory.frameNum;
        trajectoryStart = trajectoryEnd - trajectoryLength + 1; 
        
        for k = 1 : trajectoryLength
            
            currentFrame = trajectoryStart + k - 1;
            descriptor = [ currentTrajectory.trajectory(1,k); ...
                           currentTrajectory.trajectory(2,k); ...
                           currentTrajectory.HOF(:) ];
            
%             if ~isempty(trajectoriesByFrames{currentFrame}); trajectoriesByFrames{currentFrame} = []; end
            trajectoriesByFrames{currentFrame} = [trajectoriesByFrames{currentFrame} , descriptor];

        end
        
    end
    
end