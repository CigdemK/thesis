clear,clc;
warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% VARIABLE DEFINITIONS
FLAG_FIRST_NOISE = 0;

SCREEN_SIZE = [254,260];

SUBJECT_DATA = {'E:/cemre_events/elifpart1.txt',...
                'E:/cemre_events/fatihpart1.txt',...
                'E:/cemre_events/gunaypart1.txt'};

FOLDER_NAME = 'F:\Thesis\Hollywood2-actions\Hollywood2\AVIClips\';
CLIP_NAMES_PART2 = {'actioncliptest00001.avi', 'actioncliptest00004.avi', ...
                    'actioncliptest00010.avi', 'actioncliptest00021.avi', ...
                    'actioncliptest00031.avi', 'actioncliptest00057.avi', ...
                    'actioncliptest00078.avi', 'actioncliptest00084.avi', ...
                    'actioncliptest00340.avi', 'actioncliptest00351.avi', ...
                    'actioncliptest00710.avi', 'actioncliptest00734.avi', ...
                    'actioncliptest00795.avi', 'actioncliptest00882.avi', ...
                    'actioncliptrain00059.avi', 'actioncliptrain00066.avi',...
                        'actioncliptrain00763.avi_oos', 'actioncliptrain00684.avi_oos', ...
                        'actioncliptrain00239.avi_oos', 'actioncliptrain00242.avi_oos', ...
                        'actioncliptrain00300.avi_oos', 'actioncliptrain00308.avi_oos', ...
                        'actioncliptrain00438.avi_oos', 'actioncliptrain00463.avi_oos', ...
                        'actioncliptrain00484.avi_oos', 'actioncliptrain00491.avi_oos', ...
                        'actioncliptrain00546.avi_oos', 'actioncliptrain00553.avi_oos', ...
                        'actioncliptrain00598.avi_oos', 'actioncliptrain00626.avi_oos', ...
                        'actioncliptrain00629.avi_oos'};
CLIP_NAMES_PART1 = {'actioncliptrain00239.avi', 'actioncliptrain00242.avi', ...
                    'actioncliptrain00300.avi', 'actioncliptrain00308.avi', ...
                    'actioncliptrain00438.avi', 'actioncliptrain00463.avi', ...
                    'actioncliptrain00484.avi', 'actioncliptrain00491.avi', ...
                    'actioncliptrain00546.avi', 'actioncliptrain00553.avi', ...
                    'actioncliptrain00598.avi', 'actioncliptrain00626.avi', ...
                    'actioncliptrain00629.avi', 'actioncliptrain00684.avi', ...
                    'actioncliptrain00763.avi', ...
                        'actioncliptrain00066.avi_oos', ...
                        'actioncliptrain00059.avi_oos', 'actioncliptest00882.avi_oos', ...
                        'actioncliptest00795.avi_oos', 'actioncliptest00734.avi_oos', ...
                        'actioncliptest00710.avi_oos', 'actioncliptest00351.avi_oos', ...
                        'actioncliptest00340.avi_oos', 'actioncliptest00084.avi_oos', ...
                        'actioncliptest00078.avi_oos', 'actioncliptest00057.avi_oos', ...
                        'actioncliptest00031.avi_oos', 'actioncliptest00021.avi_oos', ...
                        'actioncliptest00010.avi_oos', 'actioncliptest00004.avi_oos', ...
                        'actioncliptest00001.avi_oos'};
                
%% INITIALIZATION

numberOfClips = size(CLIP_NAMES_PART1,2);
numberOfSubjects = size(SUBJECT_DATA,2);

clips = cell(numberOfClips,1);
clip_info = cell(numberOfClips,1);

for i = 1:numberOfClips
    
    clips{i} = cell(numberOfSubjects,1);
    
    for k = 1:numberOfSubjects
        clips{i}{k} = zeros(1500,3);
        clips{i}{k}(1,1) = 2; % keep the index of the next available row in the first element
    end
    
end
tic; 

%% READ VIDEOS AND SAVE FRAMES

for i = 1:numberOfClips
    
    regexres=regexp(CLIP_NAMES_PART1{i},'_','split');
    video = VideoReader([FOLDER_NAME regexres{1}]);
    clip_info{i} = [video.Height,video.Width,video.FrameRate,video.NumberOfFrames];
 
    frames = read(video);
    regexres=regexp(CLIP_NAMES_PART1{i},'.avi','split');
%     mkdir([FOLDER_NAME regexres{1}]);
    for k = 1:video.NumberOfFrames
        imwrite( frames(:,:,:,k),[[FOLDER_NAME regexres{1}] '\LeMeurOOSStats\frame_' num2str(k) '.bmp'] , 'bmp');
        imwrite( frames(:,:,:,k),[[FOLDER_NAME regexres{1}] '\LeMeurNormalStats\frame_' num2str(k) '.bmp'] , 'bmp');
    end    
    toc;

end

%% READ USER DATA AS OUTPUTTED FROM E-PRIME

for i = 1:numberOfSubjects
    
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
        if(currentSubjectData{k}(191)== '0')
            FLAG_FIRST_NOISE = 1;
            continue;
        end
        
        % when XDAT == 1 by default but event 1 is not happening
        if FLAG_FIRST_NOISE 
            currentClip = str2double(currentSubjectData{k}(191:194)); % assumption: number of clips < 1000

            nextAvailableRow = clips{currentClip}{i}(1,1);

            timeStamp = str2num(currentSubjectData{k}(71:79)); %timestamp
            clips{currentClip}{i}(nextAvailableRow,1) = timeStamp(1) + 0.001 * timeStamp(2) ;
            clips{currentClip}{i}(nextAvailableRow,2) = round(str2double(currentSubjectData{k}(231:233)) / SCREEN_SIZE(2) * clip_info{currentClip}(2)); % horizontal coords
            clips{currentClip}{i}(nextAvailableRow,3) = round(str2double(currentSubjectData{k}(251:253)) / SCREEN_SIZE(1) * clip_info{currentClip}(1)); % vertical coords

            clips{currentClip}{i}(1,1) = clips{currentClip}{i}(1,1) + 1;
        end
    end
   
    FLAG_FIRST_NOISE = 0;
end
toc;

%% SAVE DATA AS INPUTS FOR LEMEUR

for i = 1:numberOfClips
    
    nFrames = clip_info{i}(4);
    
    regexres=regexp(CLIP_NAMES_PART1{i},'.avi','split');
    
    if regexp(CLIP_NAMES_PART1{i},'oos');
        mode = 'OOS';
    else
        mode = 'Normal';
    end
    mkdir([FOLDER_NAME regexres{1} '\LeMeur' mode 'Stats']);
    
    for k = 1:nFrames

        indices = round((k-1)/clip_info{i}(3)/0.017)+2:round((k)/clip_info{i}(3)/0.017)+1;
        strToWrite = [];
        
        for subject = 1:numberOfSubjects
            newLine = '';
            for ind = indices
                x = clips{i}{subject}(ind,2);
                y = clips{i}{subject}(ind,3);
                if x > 0 && y > 0 
                    newLine = [newLine num2str(x) ' ' num2str(y) ' 17 '] ;
                end
            end
            
            if(subject ~= numberOfSubjects)
                if ~strcmp(newLine,'')
                    strToWrite = [strToWrite newLine '-1 -1 -1\n'];
                end
            else
                strToWrite = [strToWrite newLine '-1 -1 -1'];
            end
            
        end
 
        fid = fopen( [FOLDER_NAME regexres{1} '\LeMeur' mode 'Stats\frame_' num2str(k) '.stat'], 'wt' );
        fprintf(fid,strToWrite);
        fclose(fid);
        
    end
    toc;
end




