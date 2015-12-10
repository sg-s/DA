% makes a control pradigm for flcikering stim on top of a background

% some parameters
sr = 1e4; % sampling rate
b = [.1 .2 .5 .8 1]; % background odor concentration (0-5)
a_min = 0; % min foreground odor concentration (0-5)
a_max= 1.5; % max foreground odor conc
a_total = 2; % total flow from the two secondary MFCs. make sure a_total << main_air (in flow units)
             % also make sure a_total > a_max, otherwise -ve values will be
             % written to the MFC
w_min = 100;  %$ nibnimum pulse width
w_max = 3000; % maximum pulse width
main_air = 2; % main air flow rate (0-5, in units of L/min)

t_on = 5;
T = 30;
t_off = T-5;

contrast = (a_max+a_min)./(2*b); % this is the weber contrast 

% make one paradigm initially
clear ControlParadigm
ControlParadigm.Name = strcat('contrast=',oval(contrast(1),2));

% make the main and background odor outputs
ControlParadigm.Outputs(1,:) = [main_air*ones(1,t_off*sr) zeros(1,(T-t_off)*sr)];
ControlParadigm.Outputs(4,:) = [b(1)*ones(1,t_off*sr) zeros(1,(T-t_off)*sr)];
ControlParadigm.Outputs(6,:) = [zeros(1,t_on*sr) ones(1,(t_off-t_on)*sr) zeros(1,(T-t_off)*sr)];

% make the foreground odor outputs
ControlParadigm.Outputs(5,:) = zeros(1,T*sr); % valve
ControlParadigm.Outputs(2,:) = zeros(1,T*sr); % Odor MFC (Mfc1)
ControlParadigm.Outputs(3,:) = zeros(1,T*sr); % clean MFC (Mfc2)

block_start = t_on*sr;
while block_start < t_off*sr
    this_pulse_height = a_min + rand*(a_max-a_min);
    this_pulse_width = floor(w_min + rand*(w_max-w_min));
    this_blank_width = floor(w_min + rand*(w_max-w_min));
    block_end = block_start+this_pulse_width+this_blank_width;
    
    % set MFC to this value for the whole block
    ControlParadigm.Outputs(2,block_start:block_end) = this_pulse_height; 
    
    % compensate for this with the other MFC
    ControlParadigm.Outputs(3,block_start:block_end) = a_total-this_pulse_height; 
    
    % open the valve
    ControlParadigm.Outputs(5,block_start+this_blank_width:block_end) = 1;
    
    % go to the next block
    block_start = block_end+1;
end


% adjust the background in each paradigm
for i = 2:length(b)
    ControlParadigm(i).Name = strcat('contrast=',oval(contrast(i),2));
    ControlParadigm(i).Outputs = ControlParadigm(1).Outputs;
    ControlParadigm(i).Outputs(4,:) = [b(i)*ones(1,t_off*sr) zeros(1,(T-t_off)*sr)];
end


% mfc1 = a + s*randn(1,T*sr);
% mfc1 = filter(ones(1,t_mfc)/t_mfc,1,mfc1);
% mfc1(mfc1<0) = 0;
% mfc1(mfc1>5) = 5;
% mfc2 = 5-mfc1;
% mfc1(1:t_on*sr) = 0; mfc1(t_off*sr:end) = 0;
% mfc2(1:t_on*sr) = 0; mfc2(t_off*sr:end) = 0;
% ControlParadigm.Outputs(2,:) = mfc1;
% ControlParadigm.Outputs(3,:) = mfc2;

% % make the valve outputs
% 
% ov = rand(1,T*sr);
% ov = filter(ones(1,t_valve)/t_valve,1,ov);
% ov(ov<.5) = 0;
% ov(ov>0)=1;
% ov(1:t_on*sr) = 0; ov(t_off*sr:end) = 0;
% ControlParadigm.Outputs(5,:)=ov;