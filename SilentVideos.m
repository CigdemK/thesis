clear;
warning('off', 'MATLAB:MKDIR:DirectoryExists');

% variable definition
NUMBER_OF_CLIPS = 31;
NUMBER_OF_SUBJECTS = 4;
SUBJECT_DATA = {'C:/Users/ckocb_000/Desktop/Cigdem/abdullah.txt', ...
              'C:/Users/ckocb_000/Desktop/Cigdem/cagri.txt', ...
              'C:/Users/ckocb_000/Desktop/Cigdem/cemre.txt', ...
              'C:/Users/ckocb_000/Desktop/Cigdem/sevgi.txt' };
CLIP_FPS=[24;24;24;24;24;24;24;24;24;24; ...
          24;24;24;24;24;24;24;24;24;24; ...
          24;24;24;24;24;24;24;24;24;24; ...
          12];
FLAG_FIRST_NOISE = 0;
LAST_CLIP_DURATION = 49;
VIDEO_SIZE = [224,548];
SCREEN_SIZE = [254,260];

% initialization
clips = cell(NUMBER_OF_CLIPS,1);
for i = 1:NUMBER_OF_CLIPS
    clips{i} = cell(NUMBER_OF_SUBJECTS,1);
    for k = 1:NUMBER_OF_SUBJECTS
        clips{i}{k} = zeros(1500,3);
        clips{i}{k}(1,1) = 2; % keep the index of the next available row in the first element
    end
end
tic; 
% read user data as outputted by e-prime
for i = 1:NUMBER_OF_SUBJECTS
    
    fid = fopen(SUBJECT_DATA{i});
    currentSubjectData = textscan(fid,'%s','Delimiter','\n');
    currentSubjectData = currentSubjectData{1};
    fclose(fid);
    
    for k = 2:size(currentSubjectData,1)
        
        % break if we reach end of data
        if strcmp(currentSubjectData{k},'') 
            break; 
        end
        
        % skip if XDAT is 0
        if(currentSubjectData{k}(181)== '0')
            FLAG_FIRST_NOISE = 1;
            continue;
        end
        
        % when XDAT == 1 by default but event 1 is not happening
        if FLAG_FIRST_NOISE 
            currentClip = str2double(currentSubjectData{k}(181:183)); % assumption: number of clips < 1000

            nextAvailableRow = clips{currentClip}{i}(1,1);

            timeStamp = str2num(currentSubjectData{k}(61:68)); %timestamp
            clips{currentClip}{i}(nextAvailableRow,1) = timeStamp(1) + 0.001 * timeStamp(2) ;
            clips{currentClip}{i}(nextAvailableRow,2) = round(str2double(currentSubjectData{k}(221:223)) / SCREEN_SIZE(2) * VIDEO_SIZE(2)); % horizontal coords
            clips{currentClip}{i}(nextAvailableRow,3) = round(str2double(currentSubjectData{k}(241:243)) / SCREEN_SIZE(1) * VIDEO_SIZE(1)); % vertical coords

            clips{currentClip}{i}(1,1) = clips{currentClip}{i}(1,1) + 1;
        end
    end
   
    FLAG_FIRST_NOISE = 0;
end
toc;
for i = 1:NUMBER_OF_SUBJECTS
    clips{end}{i}(LAST_CLIP_DURATION*61:end,:) = [];
end
toc;  
% prepare for LeMeur
for i = 1:NUMBER_OF_CLIPS
    
    lastElement = clips{i}{2}(1,1)-1;
    nFrames = ceil((clips{i}{2}(lastElement,1) - clips{i}{2}(2,1)) * CLIP_FPS(i));
    linesPerFrame = round((clips{i}{2}(1,1) - 2) / nFrames);
    
    for k = 1:nFrames
        
%         indices = round((k-1)/CLIP_FPS(i)/0.017)+2:round((k)/CLIP_FPS(i)/0.017)+1;
        indices = (k-1)*linesPerFrame+2:k*linesPerFrame+1;
        strToWrite = [];
        
        for subject = 1:NUMBER_OF_SUBJECTS
            for ind = indices
                x = clips{i}{subject}(ind,2);
                y = clips{i}{subject}(ind,3);
                strToWrite = [strToWrite num2str(x) ' ' num2str(y) ' 17 '] ;
            end
            
            if(subject ~= NUMBER_OF_SUBJECTS)
                strToWrite = [strToWrite '-1 -1 -1\n'];
            else
                strToWrite = [strToWrite '-1 -1 -1'];
            end
            
        end
        
        mkdir(['./data/clip' num2str(i)]);
        fid = fopen( ['./data/clip' num2str(i) '/frame_' num2str(k) '.stat'], 'wt' );
        fprintf(fid,strToWrite);
        fclose(fid);
        
    end
    toc;
end


