timerTic=4; % how often the timer checks

timerHandle = timer();
timerHandle.startDelay = timerTic;
timerHandle.Period = timerTic;
timerHandle.ExecutionMode = 'fixedRate';
timerHandle.TasksToExecute = inf;
timerHandle.TimerFcn = @mycallbackfunction;

start(timerHandle);