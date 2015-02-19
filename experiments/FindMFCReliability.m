% FindMFCReliability.m
% this script uses kontroller to find how reproducible the MFC is when
% driven with some exponentiated gaussian noise input, as a function of
% switching time. 

switching_time = [.05 .075 .1 .125 .2]; % in seconds
results(1).r = [];
results(1).data = [];

for i = 1:length(switching_time)
    tc_string = strcat('tc=',mat2str(switching_time(i)),';');
    ControlParadigm = MakeLargeVarianceFlicker('T=10;','s=1;',tc_string);
    
    % run it
    results(i).data = Kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',2*ones(1,10),'w',1e4);
end
