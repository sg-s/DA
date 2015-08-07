Total_Flow = 120; % mL/min, odoour flow
MFC_max = 500;
T = 30;
pulse_on = 4;
pulse_off = 5;
odor_vol = [0.2 1 2 5 10 20 30 50 120];

% make a start control paradigm
c = 1;
clear ControlParadigm
ControlParadigm(1).Name = 'start';
ControlParadigm(1).Outputs = zeros(5,10000);
% turn main air on
ControlParadigm(1).Outputs(5,:) = 1; 
c=c+1;

% make the pulse of various dilutions 
for i = 1:length(odor_vol)
    ControlParadigm(c).Name = [mat2str(odor_vol(i)),':',mat2str(Total_Flow-odor_vol(i))];
    
    ControlParadigm(c).Outputs = zeros(5,T*1e4);
    
    % set the MFCs
    ControlParadigm(c).Outputs(2,:) = (Total_Flow-odor_vol(i))/MFC_max*5;
    ControlParadigm(c).Outputs(3,:) = odor_vol(i)/MFC_max*5;
    
    % turn on the main air
    ControlParadigm(c).Outputs(5,:) = 1;
    
    % turn the first valve on
    ControlParadigm(c).Outputs(4,pulse_on*1e4:pulse_off*1e4) = 1;
    c = c+1;
    
    
end

% turn everything off
ControlParadigm(c).Name = 'end';
ControlParadigm(c).Outputs = zeros(5,10000);



save('Davide_Pulses_30sec_05_26_15_Kontroller_paradigm.mat','ControlParadigm')

