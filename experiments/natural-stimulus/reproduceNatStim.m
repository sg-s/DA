%% builds control signals to reproduce natural stimulus

%% build flash response

% parameters
flash_duration = 500; % in correct time units, 50ms
min_flow_rate = 1/10; % in V
max_flow_rate = 3; 
flow_rates = logspace(log10(min_flow_rate),log10(max_flow_rate),20);
T = 3e4;
clear ControlParadigm
for i = 1:length(flow_rates)
	ControlParadigm(i).Name = oval(flow_rates(i));
	ControlParadigm(i).Outputs = zeros(3,T);
	ControlParadigm(i).Outputs(3,:) = 1;
	ControlParadigm(i).Outputs(1,:) = flow_rates(i);
	ControlParadigm(i).Outputs(1,end) = 0;
	ControlParadigm(i).Outputs(2,2e4:2e4+500) = 1; % valve
end

disp('Building flash response curve...')

data = kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',1:length(flow_rates),'w',1e4);

peak_stim = 0*flow_rates;

for i = 1:length(data)
	peak_stim(i) = max(data(i).PID) - mean(data(i).PID(1:2e4));
end

peak_stim = peak_stim(:);
flow_rates = flow_rates(:);

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(flow_rates,peak_stim,'k+')
ff = fit(flow_rates(2:end),peak_stim(2:end),'smoothingspline');
plot(flow_rates,ff(flow_rates),'r')
xlabel('Flow Rates (V)')
ylabel('Peak PID for 50ms pulse')
set(gca,'XScale','log','YScale','log')

prettyFig('FixLogX',true);


%% Now we construct the ansatz solution for the flow rates
load('nat_stim_target.mat')

clear ControlParadigm
ControlParadigm.Name = 'Natural Stim';
ControlParadigm.Outputs = zeros(3,70e4);

% first set the air to be on always;
ControlParadigm.Outputs(3,:) = 1;

% now get the valve to turn on when we want to

for i = 1:length(ons)
	ControlParadigm.Outputs(2,ons(i):offs(i)) = 1;
end

% now, starting from the last whiff, set the MFC to the desired setpoint for 1 second
scan_x = 0:1e-3:5;
for i = length(ons):-1:1
	[~,x] = min(abs(ff(scan_x)-peak_S(i)));
	ControlParadigm.Outputs(1,ons(i)-1e4:ons(i)) = scan_x(x); 
end

save('Ansatz_Nat_Stim_kontroller_paradigm.mat','ControlParadigm')





