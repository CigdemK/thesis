function mycallbackfunction(~,~)
    if feature('IsDebugMode') % undocumented thanks to CatzLoveJazz
        load handel
        sound(y,Fs);
        evalin('base','stop(timerHandle)') % stop the timer
    end
end


% stop(timerHandle)
% delete(timerHandle)

playing = false
function mycallbackfunction(~,~)
    if feature('IsDebugMode')
    	if !playing % undocumented thanks to CatzLoveJazz
        	load handel
        	sound(y,Fs);
        	playing = true
	end	
    else
    	playing = false
   	end

end
