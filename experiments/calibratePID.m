%% calibratePID.m
% calibrates PID using the depletion technique 

% choose flow rate
MFC_setpoint = .75; % V

% make control paradigms
ControlParadigm.Name = 'Zero';
ControlParadigm.Outputs = zeros(3,10);
ControlParadigm.Outputs(3,:) = 1;


% measure the zero point
first_blank = kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',ones(5,1),'w',10);

PID = [];
t = [];
this_PID = Inf;
i = 1;

% make the figure and the plots
f = figure; hold on
h = plot(NaN,NaN,'k+');

% make control paradigms
ControlParadigm.Name = 'Flow';
ControlParadigm.Outputs = ones(3,10);
ControlParadigm.Outputs(1,:) = MFC_setpoint;

while this_PID > mean(mean([first_blank.PID]')) + 3*std(mean([first_blank.PID]'))
    clc
    data = kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',1,'w',10);
    this_PID = mean(data.PID(:));
    PID = [PID; this_PID];
    t = [t; now];
    set(h,'XData',t,'YData',PID);
    datetick('x')
end