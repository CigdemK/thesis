function MyCallbackFunction(~,~)
    
    if ~exits(flagPlaying)
        global flagPlaying = false;
    end

    if feature('IsDebugMode')
    	if ~flagPalying
        	load handel;
        	sound(y,Fs);
        	flagPlaying = true;
        end
    else
    	flagPlaying = false;
    end
    
end

% 
% function mycallbackfunction(~,~)
%     if feature('IsDebugMode') % undocumented thanks to CatzLoveJazz
%         load handel
%         sound(y,Fs);
%         evalin('base','stop(timerHandle)') % stop the timer
%     end
% end
% 
% 
% % stop(timerHandle)
% % delete(timerHandle)


