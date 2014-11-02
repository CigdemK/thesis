function EnableDebugModeNotification(~,~)
    
    if ~exist('flagPlaying', 'var')
        evalin('base','global flagPlaying');
        flagPlaying = false;
    end

    if feature('IsDebugMode')
    	if ~flagPlaying
        	load handel;
        	sound(y,Fs);
        	flagPlaying = true;
        end
    else
    	flagPlaying = false;
    end
    
end