clear;
clc;
Total_Flow = 120; % mL/min, odoour flow
MFC_max = 500;
T = 30;
pulse_on = 4;
pulse_off = 5;
% odor_vol = [0.2 1 2 5 10 20 30 50 120];
odor_vol = [1 10 50 120];

% make a start control paradigm
c = 1;
clear ControlParadigm
ControlParadigm(1).Name = 'start';
% 1dilute, 2light, 3mfc500, 4valve, 5plit, 6main
ControlParadigm(1).Outputs = zeros(6,10000); 
% turn main air on
ControlParadigm(1).Outputs(6,:) = 1; 
c=c+1;

% make split switches to check air balance
split_width = 2; % sec
ControlParadigm(c).Name = ['main_air_split:',mat2str(split_width),'_sec'];

ControlParadigm(c).Outputs = zeros(6,T*1e4);

% turn on the main air
ControlParadigm(c).Outputs(6,:) = 1;


 % pulse the split valve

 for k = split_width:2*split_width:int64(T);
     if (k+split_width)<T
        ControlParadigm(c).Outputs(5,k*1e4:(k+split_width)*1e4) = 1;
     end
 end
c = c+1;


% make the pulse of various dilutions without air split 
for i = 1:length(odor_vol)
    ControlParadigm(c).Name = [mat2str(odor_vol(i)),':',mat2str(Total_Flow-odor_vol(i))];
    
    ControlParadigm(c).Outputs = zeros(6,T*1e4);
    
    % set the MFCs
    ControlParadigm(c).Outputs(1,:) = (Total_Flow-odor_vol(i))/MFC_max*5;
    ControlParadigm(c).Outputs(3,:) = odor_vol(i)/MFC_max*5;
    
    % turn on the main air
    ControlParadigm(c).Outputs(6,:) = 1;
    
    % turn the first valve on
    ControlParadigm(c).Outputs(4,pulse_on*1e4:pulse_off*1e4) = 1;

    c = c+1;
end

% make the pulse of various dilutions 
for i = 1:length(odor_vol)
    ControlParadigm(c).Name = ['air_split_', mat2str(odor_vol(i)),':',mat2str(Total_Flow-odor_vol(i))];
    
    ControlParadigm(c).Outputs = zeros(6,T*1e4);
    
    % set the MFCs
    ControlParadigm(c).Outputs(1,:) = (Total_Flow-odor_vol(i))/MFC_max*5;
    ControlParadigm(c).Outputs(3,:) = odor_vol(i)/MFC_max*5;
    
    % turn on the main air
    ControlParadigm(c).Outputs(6,:) = 1;
    
    % turn the first valve on
    ControlParadigm(c).Outputs(4,pulse_on*1e4:pulse_off*1e4) = 1;
    
    
     % turn the split valve after the valve off
    ControlParadigm(c).Outputs(5,pulse_off*1e4:end-1) = 1;
    c = c+1;
end


% turn everything off
ControlParadigm(c).Name = 'end';
ControlParadigm(c).Outputs = zeros(6,10000);
ControlParadigm(c-1).Outputs(3,end) = 0;



save('Davide_Pulses_30sec_split_06_01_15_Kontroller_paradigm.mat','ControlParadigm')

