%% analysis of random stim

clear all
filename = '/data/random-stim/final_2011_06_14_ab3A_1o3ol3X-3_20ml_30sec_30ms_rand 2.mat'
deltat = 10^(-4);

num_files = 1;
t_ini = 9;
t_f = 35;
win = 0.03;
sliding = 0.003;
memory = round(1/sliding);
shift_input = round(0.1/sliding);
    
for jj = 1:num_files
    load(filename);
    trialsPID = size(PID,2);
    trials = size(stA,2);
    
    % bin and average
    [timesr PIDs sr msr sdr] = bin_traces(PID, stA, deltat, win, sliding);
    time = timesr(find(timesr>t_ini,1):find(timesr>t_f,1));   
    
    % PID
    meanPID = mean(squeeze(PIDs(1,:,:)), 1);      
    sdPID = std(squeeze(PIDs(1,:,:)), 1);      
%     
%     ptemp = meanPID(find(timesr>t_ini,1):find(timesr>t_f,1));
%     PID_offset(:,jj) = mean(ptemp);
%     ptemp =  ptemp - smooth( ptemp, 5/sliding)';
%     ptemp = ptemp + PID_offset(:,jj);
%     meanPID(find(timesr>t_ini,1):find(timesr>t_f,1)) =  ptemp;
%     figure; hold on
%     plot(meanPID, 'k')
%     meanPID = butter_filter(meanPID, .001, 0, sliding, 3, 3);
%     plot(meanPID, 'b')

%     rasterplot(stim_signal, PIDs, meanPID, sdPID/sqrt(trialsPID),  stA, msr, sdr/sqrt(trials), deltat, timesr, 25, 35)
    p(:,jj) = meanPID(find(timesr>t_ini,1):find(timesr>t_f,1)); 
%     p(:,jj) = detrend(p(:,jj));
    psd(:,jj) = sdPID(find(timesr>t_ini,1):find(timesr>t_f,1));
%     if jj==2
%         p(:,jj) = p(:,jj)*.4 + 0.004;
%     end      
    pmax(:,jj) = max(p(:,jj));
    p(:,jj) = p(:,jj)/pmax(:,jj);
    psd(:,jj) = psd(:,jj)/pmax(:,jj); 
    
    % ORN
    f(:,jj) = msr(find(timesr>t_ini,1):find(timesr>t_f,1));%squeeze(sr(1,1,find(timesr>t_ini,1):find(timesr>t_f,1))); % 
    fsd(:,jj) = sdr(find(timesr>t_ini,1):find(timesr>t_f,1));%zeros(length(f(:,jj)),1); % 
    
    % valve
    stim_signal = round(stim_signal);
    [timesr stimtemp] = sliding_average(squeeze(stim_signal(1,1,:)), deltat, win, sliding);
    stim(:,jj) = stimtemp(find(timesr>t_ini,1):find(timesr>t_f,1));    
    stim(:,jj) = stim(:,jj) - mean(stim(:,jj));
    stim(:,jj) = stim(:,jj)/max(stim(:,jj));
   
    % LN model for PID response from valve
    fitfun = 0; % type of non-linear static function: 0==(linear fit;) 1==(hill function)
    figon = 1; % 1==(shows plot of LN model and crossprediction)
    reg = 500; %[0 10.^[-1:0.5:3]];
    regL2 = 1; % type of regularization: 0==(elastic net); 1==(L2 regularization by me)
    [input_s2p(:,jj), f2_s2p(:,jj), fSTA_s2p(:,jj),  fLN_s2p(:,jj), a_s2p(:,jj), D_s2p(:,jj), NR_s2p(:,jj), regopt_s2p(:,jj), lambda_s2p(:,jj)] = build_LNmodel_optreg(time, stim(:,jj), p(:,jj), psd(:,jj), trialsPID, memory, shift_input, sliding, reg, fitfun, figon, regL2);

    % LN model for ORN response from valve
    fitfun = 1; % type of non-linear static function: 0==(linear fit;) 1==(hill function)
    figon = 1; % 1==(shows plot of LN model and crossprediction)
    reg = 1000; %[10.^[-1:0.5:3]];
    regL2 = 1; % type of regularization: 0==(elastic net); 1==(L2 regularization by me)
    [input_s2n(:,jj), f2_s2n(:,jj), fSTA_s2n(:,jj), fLN_s2n(:,jj), a_s2n(:,jj), D_s2n(:,jj), NR_s2n(:,jj), regopt_s2n(:,jj), lambda_s2n(:,jj)]= build_LNmodel_optreg(time, stim(:,jj), f(:,jj), fsd(:,jj), trials, memory, shift_input, sliding, reg, fitfun, figon, regL2);
  
    % LN model for ORN response from PID
    fitfun = 1; % type of non-linear static function: 0==(linear fit;) 1==(hill function)
    figon = 1; % 1==(shows plot of LN model and crossprediction)
    reg = 50; %[10.^[-1:0.5:3]];
    regL2 = 1; % type of regularization: 0==(elastic net); 1==(L2 regularization by me)
    [input_p2n(:,jj), f2_p2n(:,jj), fSTA_p2n(:,jj), fLN_p2n(:,jj), a_p2n(:,jj), D_p2n(:,jj), NR_p2n(:,jj), regopt_p2n(:,jj), lambda_p2n(:,jj)]= build_LNmodel_optreg(time, p(:,jj), f(:,jj), fsd(:,jj), trials, memory, shift_input, sliding, reg, fitfun, figon, regL2);

  
end
%% plots
color = [1 0 0; 0 .75 .75; 0 .75 .75];
t_shifted = (1:memory+1)*sliding - 0.1;

figure; hold on
subplot(3,3,1); hold on
for jj = 1:num_files
    plot(t_shifted, a_s2p(:,jj), 'linewidth', 2, 'color', color(jj,:))
end
xlim([t_shifted(1) t_shifted(end)])

subplot(3,3,4); hold on
for jj = 1:num_files
    plot(t_shifted, a_s2n(:,jj), 'linewidth', 2, 'color', color(jj,:))
end
xlim([t_shifted(1) t_shifted(end)])
ylim([-1 1.2])

subplot(3,3,7); hold on
for jj = 1:num_files
    plot(t_shifted+0.003, a_p2n(:,jj), 'linewidth', 2, 'color', color(jj,:))
end
xlim([t_shifted(1) t_shifted(end)])
ylim([-1 1.2])

% non-linear part
Ltot = length(f2_s2p(:,jj));
Ltest = floor(Ltot/3);
itest = (Ltot-Ltest+1:Ltot)';

for jj=1:num_files
    subplot(3,3,1+jj); hold on
    plot(fSTA_s2p(itest,jj), f2_s2p(itest,jj), '.', 'color', [.8, .8, .8])
    D = D_s2p(:,jj);
    ffitted = @(x) D(1)*x + D(2);
    fplot(ffitted, [min(fSTA_s2p(itest,jj)) max(fSTA_s2p(itest,jj))]);
end

for jj=1:num_files
    subplot(3,3,4+jj); hold on
    plot(fSTA_s2n(itest,jj), f2_s2n(itest,jj), '.', 'color', [.8, .8, .8])
    D = D_s2n(:,jj);
    ffitted = @(x) D(1)*(1./(1+(D(2)./(x-D(4))).^D(3)));
    fplot(ffitted, [min(fSTA_s2n(itest,jj)) max(fSTA_s2n(itest,jj))]);
end

for jj=1:num_files
    subplot(3,3,7+jj); hold on
    plot(fSTA_p2n(itest,jj), f2_p2n(itest,jj), '.', 'color', [.8, .8, .8])
    D = D_p2n(:,jj);
    ffitted = @(x) D(1)*(1./(1+(D(2)./(x-D(4))).^D(3)));
    fplot(ffitted, [min(fSTA_p2n(itest,jj)) max(fSTA_p2n(itest,jj))]);
end


%% fitting kernel p2n with laguerre polynomials

t_ini_all = [32 32 34]; % guess the delay 

t_ini_all = t_ini_all-1;

D_temp = [];
Da = [];
for jj=1:num_files
    
    t_ini = t_ini_all(jj);
    delay = (t_ini)*sliding;

    tt = (t_ini:memory+1)'*sliding - delay;
    aa = a_p2n(t_ini:end,jj);
    % tt(aa==0)=[];
    % aa(aa==0)=[];
    laguerre2 = fittype('exp(-l1*t)*b1*(cos(l4*t)+sin(l4*t)) + exp(-l2*t)*b2*(cos(l3*t)+sin(l3*t))', 'coefficients', {'l1', 'b1', 'l2', 'b2', 'l3', 'l4'}, 'independent', 't');

    options = fitoptions('Method', 'NonlinearLeastSquares');
    options.MaxFunEval = 5000;
    options.MaxIter = 1000;
%     options.StartPoint = [300 -10 50 -10 -20 50]; % both
    if jj==1||2
        options.StartPoint = [500 -10 10 10 .01 10]; % 1but
    elseif jj==3
        options.StartPoint = [100 -10 10 10 -10 -1]; % 1o3ol
    end
    fitN = fit(tt, aa, laguerre2, options);
    Da(:,jj) = coeffvalues(fitN);
    
%     fa = @(t) exp(-Da(1)*(t-delay))*Da(2)*(cos(Da(6)*(t-delay))+sin(Da(6)*(t-delay))) + exp(-Da(3)*(t-delay))*Da(4)*(cos(Da(5)*(t-delay))+sin(Da(5)*(t-delay)));
%     figure
%     fplot(fa, [0 .6])
%     t0 = fzero(fa,delay);
    a_temp(:,jj) = [zeros(t_ini-1,1); exp(-Da(1).*tt)*Da(2).*(cos(Da(6)*tt)+sin(Da(6)*tt)) + exp(-Da(3).*tt)*Da(4).*(cos(Da(5)*tt)+sin(Da(5)*tt))];

    tt = [(1:t_ini-1)'*sliding; tt+delay];
    
    stim_temp = p(shift_input:end, jj);
    L = length(stim_temp);
    s = zeros(L-memory, memory+1);
    for i=1:L-memory
        s(i,:) = stim_temp(memory+i:-1:i);
    end
    fSTA_temp = s*a_temp(:,jj);
    
    figure
    subplot(2,2,3); hold on
    plot(tt, a_p2n(:,jj), '.k')
    plot(tt, a_temp(:,jj), '-r', 'linewidth', 2)
    xlim([0 .6])
    subplot(2,2,4); hold on
    plot(fSTA_temp, f(memory+1:end-shift_input+1,jj), '.', 'color', [.8 .8 .8])
    
    fun_fit = fittype('A*(1/(1+(K/x)^n))', 'coefficients',{'A' 'K', 'n'}, 'independent', 'x');
    options = fitoptions('Method', 'NonlinearLeastSquares');
    options.StartPoint = [100 mean(fSTA_temp-min(fSTA_temp)) 3];
    fitN = fit(fSTA_temp-min(fSTA_temp), f(memory+1:end-shift_input+1,jj), fun_fit, options);
    D_temp(jj,:) = [coeffvalues(fitN) min(fSTA_temp)];
    ffitted = @(x) D_temp(jj,1)*(1/(1+(D_temp(jj,2)/(x-D_temp(jj,4)))^D_temp(jj,3)));
    fplot(ffitted, [min(fSTA_temp) max(fSTA_temp)], 'r');
    fLN_temp = D_temp(jj,1)*(1./(1+(D_temp(jj,2)./(fSTA_temp-D_temp(jj,4))).^D_temp(jj,3)));
    
    subplot(2,2,1:2); hold on    
    plot(f(memory+1:end-shift_input+1,jj), 'k', 'linewidth', 1)
    plot(fLN_temp, 'r', 'linewidth', 1)
    plot(fLN_p2n(:,jj), 'g', 'linewidth', 1)
    xlim([0 3000])
    
    Pr = mean(abs(f(memory+1:end-shift_input+1,jj)-fLN_temp).^2);
    Pn = mean(fsd(memory+1:end-shift_input+1,jj).^2);
    NR_temp(jj) = sqrt(Pr/Pn);
end

%%
figure
for jj = 1:num_files
    subplot(2,3,jj); hold on
    plot(tt, a_p2n(:,jj), '.k')
    plot(tt, a_temp(:,jj), '-r', 'linewidth', 2)
    xlim([0 .6])
    ylim([-.6 1.1])
end
subplot(2,3,4)
bar([NR_p2n' NR_temp'])

subplot(2,3,5)
plot(tt, a_temp, 'linewidth', 2)
xlim([0 .6])
ylim([-.6 1.1])

subplot(2,3,6)
plot(tt, a_p2n, 'linewidth', 2)
xlim([0 .6])
ylim([-.6 1.1])

%% fitting kernel s2n with laguerre polynomials

t_ini_all = [50 52]; % guess the delay 
t_ini_all = t_ini_all-1;

D_temp = [];

for jj=1:num_files
    
    t_ini = t_ini_all(jj);
    delay = (t_ini)*sliding;

    tt = (t_ini:memory+1)'*sliding - delay;
    aa = a_s2n(t_ini:end,jj);
    % tt(aa==0)=[];
    % aa(aa==0)=[];
    laguerre2 = fittype('exp(-l1*t)*b1*(cos(l4*t)+sin(l4*t)) + exp(-l2*t)*b2*(cos(l3*t)+sin(l3*t))', 'coefficients', {'l1', 'b1', 'l2', 'b2', 'l3', 'l4'}, 'independent', 't');
    options = fitoptions('Method', 'NonlinearLeastSquares');
    options.MaxFunEval = 1000;
    options.StartPoint = [200 -10 10 10 -.01 10]; %[20 0 100 -1000];
    fitN = fit(tt, aa, laguerre2, options);
    Da = coeffvalues(fitN);
%     figure; hold on
%     plot(tt, aa, 'ok')
%     ffitted = @(t) exp(-Da(1)*t)*Da(2) + exp(-Da(3)*t)*(Da(4)*cos(Da(6)*t)+Da(5)*sin(Da(6)*t));
%     fplot(ffitted, [0 (memory+1)*sliding], 'r');
    
    a_temp(:,jj) = [zeros(t_ini-1,1); exp(-Da(1).*tt)*Da(2).*(cos(Da(6)*tt)+sin(Da(6)*tt)) + exp(-Da(3).*tt)*Da(4).*(cos(Da(5)*tt)+sin(Da(5)*tt))];
    tt = [(1:t_ini-1)'*sliding; tt+delay];
    
    stim_temp = stim(shift_input:end, jj);
    L = length(stim_temp);
    s = zeros(L-memory, memory+1);
    for i=1:L-memory
        s(i,:) = stim_temp(memory+i:-1:i);
    end
    fSTA_temp = s*a_temp(:,jj);
    
    figure
    subplot(2,2,3); hold on
    plot(tt, a_s2n(:,jj), '.k')
    plot(tt, a_temp(:,jj), '-r', 'linewidth', 2)
    xlim([0 .6])
    subplot(2,2,4); hold on
    plot(fSTA_temp, f(memory+1:end-shift_input+1,jj), '.', 'color', [.8 .8 .8])
    
    fun_fit = fittype('A*(1/(1+(K/x)^n))', 'coefficients',{'A' 'K', 'n'}, 'independent', 'x');
    options = fitoptions('Method', 'NonlinearLeastSquares');
    options.StartPoint = [100 mean(fSTA_temp-min(fSTA_temp)) 3];
    fitN = fit(fSTA_temp-min(fSTA_temp), f(memory+1:end-shift_input+1,jj), fun_fit, options);
    D_temp(jj,:) = [coeffvalues(fitN) min(fSTA_temp)];
    ffitted = @(x) D_temp(jj,1)*(1/(1+(D_temp(jj,2)/(x-D_temp(jj,4)))^D_temp(jj,3)));
    fplot(ffitted, [min(fSTA_temp) max(fSTA_temp)], 'r');
    fLN_temp = D_temp(jj,1)*(1./(1+(D_temp(jj,2)./(fSTA_temp-D_temp(jj,4))).^D_temp(jj,3)));
    
    subplot(2,2,1:2); hold on    
    plot(f(memory+1:end-shift_input+1,jj), 'k', 'linewidth', 1)
    plot(fLN_temp, 'r', 'linewidth', 1)
    plot(fLN_s2n(:,jj), 'g', 'linewidth', 1)
    xlim([0 3000])
    
    Pr = mean(abs(f(memory+1:end-shift_input+1,jj)-fLN_temp).^2);
    Pn = mean(fsd(memory+1:end-shift_input+1,jj).^2);
    NR_temp(jj) = sqrt(Pr/Pn);
end
%%
figure
for jj = 1:num_files
    subplot(2,3,jj); hold on
    plot(tt, a_s2n(:,jj), '.k')
    plot(tt, a_temp(:,jj), '-r', 'linewidth', 2)
    xlim([0 .6])
    ylim([-.6 1.1])
end
subplot(2,3,4)
bar([NR_s2n' NR_temp'])

subplot(2,3,5)
plot(tt, a_temp, 'linewidth', 2)
xlim([0 .6])
ylim([-.6 1.1])

subplot(2,3,6)
plot(tt, a_s2n, 'linewidth', 2)
xlim([0 .6])
ylim([-.6 1.1])

%% fitting kernel s2p with laguerre polynomials

t_ini_all = [47 47]; % guess the delay 30=p2n 49/50=s2n 48/45=s2p
t_ini_all = t_ini_all-1;

D_temp = [];

for jj=1:num_files
    
    t_ini = t_ini_all(jj);
    delay = (t_ini)*sliding;

    tt = (t_ini:memory+1)'*sliding - delay;
    aa = a_s2p(t_ini:end,jj);
    % tt(aa==0)=[];
    % aa(aa==0)=[];
    laguerre2 = fittype('exp(-l1*t)*b1 + exp(-l2*t)*b2 + exp(-l3*t)*b3', 'coefficients', {'l1', 'b1', 'l2', 'b2', 'l3', 'b3'}, 'independent', 't');
    options = fitoptions('Method', 'NonlinearLeastSquares');
    options.MaxFunEval = 1000;
    options.StartPoint = [200 -10 30 -10 1 -10]; %[20 0 100 -1000];
    fitN = fit(tt, aa, laguerre2, options);
    Da = coeffvalues(fitN);
%     figure; hold on
%     plot(tt, aa, 'ok')
%     ffitted = @(t) exp(-Da(1)*t)*Da(2) + exp(-Da(3)*t)*(Da(4)*cos(Da(6)*t)+Da(5)*sin(Da(6)*t));
%     fplot(ffitted, [0 (memory+1)*sliding], 'r');
    
    a_temp(:,jj) = [zeros(t_ini-1,1); exp(-Da(1).*tt)*Da(2) + exp(-Da(3).*tt)*Da(4) + exp(-Da(5).*tt)*Da(6)];
    tt = [(1:t_ini-1)'*sliding; tt+delay];
    
    stim_temp = stim(shift_input:end, jj);
    L = length(stim_temp);
    s = zeros(L-memory, memory+1);
    for i=1:L-memory
        s(i,:) = stim_temp(memory+i:-1:i);
    end
    fSTA_temp = s*a_temp(:,jj);
    
    figure
    subplot(2,2,3); hold on
    plot(tt, a_s2p(:,jj), '.k')
    plot(tt, a_temp(:,jj), '-r', 'linewidth', 2)
    xlim([0 1])
    subplot(2,2,4); hold on
    plot(fSTA_temp, p(memory+1:end-shift_input+1,jj), '.', 'color', [.8 .8 .8])
      
    fun_fit = fittype('poly1');
    fitN = fit(fSTA_temp, p(memory+1:end-shift_input+1,jj), fun_fit);
    D = coeffvalues(fitN);
    ffitted = @(x) D(1)*x + D(2);
    fplot(ffitted, [min(fSTA_temp) max(fSTA_temp)], 'r');
    fLN_temp = D(1)*fSTA_temp + D(2);
    
    subplot(2,2,1:2); hold on    
    plot(p(memory+1:end-shift_input+1,jj), 'k', 'linewidth', 1)
    plot(fLN_temp, 'r', 'linewidth', 1)
    plot(fLN_s2p(:,jj), 'g', 'linewidth', 1)
    xlim([0 3000])
    
    Pr = mean(abs(p(memory+1:end-shift_input+1,jj)-fLN_temp).^2);
    Pn = mean(psd(memory+1:end-shift_input+1,jj).^2);
    NR_temp(jj) = sqrt(Pr/Pn);
end
%%
figure
for jj = 1:num_files
    subplot(2,3,jj); hold on
    plot(tt, a_s2p(:,jj), '.k')
    plot(tt, a_temp(:,jj), '-r', 'linewidth', 2)
    xlim([0 1])
    ylim([-.6 1.1])
end
subplot(2,3,4)
bar([NR_s2p' NR_temp'])

subplot(2,3,5)
plot(tt, a_temp, 'linewidth', 2)
xlim([0 1])
ylim([-.6 1.1])

subplot(2,3,6)
plot(tt, a_s2p, 'linewidth', 2)
xlim([0 1])
ylim([-.6 1.1])


%% cross-prediction

clear all
files = dir('final_*');
deltat = 10^(-4);

files = files([2 3]); %[5 3 6 1 2 4]);
num_files = size(files,1);
t_ini = 9;
t_f = 35;
win = 0.03;
sliding = 0.003;
memory = round(1/sliding);
shift_input = round(0.1/sliding);
    
for jj = 1:num_files
    load(files(jj).name); display(files(jj).name) 
    trialsPID = size(PID,2);
    trials = size(stA,2);
    
    % bin and average
    [timesr PIDs sr msr sdr] = bin_traces(PID, stA, deltat, win, sliding);
    time = timesr(find(timesr>t_ini,1):find(timesr>t_f,1));   
    
    % PID
    meanPID = mean(squeeze(PIDs(1,:,:)), 1);      
    sdPID = std(squeeze(PIDs(1,:,:)), 1); 
    
%     ptemp = meanPID(find(timesr>t_ini,1):find(timesr>t_f,1));
%     PID_offset(:,jj) = mean(ptemp);
%     ptemp =  ptemp - smooth( ptemp, 5/sliding)';
%     ptemp = ptemp + PID_offset(:,jj);
%     meanPID(find(timesr>t_ini,1):find(timesr>t_f,1)) =  ptemp;

%     rasterplot(stim_signal, PIDs, meanPID, sdPID/sqrt(trialsPID),  stA, msr, sdr/sqrt(trials), deltat, timesr, 20, 30)
    p(:,jj) = meanPID(find(timesr>t_ini,1):find(timesr>t_f,1)); 
    PID_offset(:,jj) = mean(p(:,jj));
    p(:,jj) = detrend(p(:,jj))+PID_offset(:,jj);
    psd(:,jj) = sdPID(find(timesr>t_ini,1):find(timesr>t_f,1));
%     if jj==2
%         p(:,jj) = p(:,jj)*.4 + 0.004;
%     end      
%     pmax(:,jj) = max(p(:,jj));
%     p(:,jj) = p(:,jj)/pmax(:,jj);
%     psd(:,jj) = psd(:,jj)/pmax(:,jj); 
    
    % valve
    stim_signal = round(stim_signal);
    [timesr stimtemp] = sliding_average(squeeze(stim_signal(1,1,:)), deltat, win, sliding);
    stim(:,jj) = stimtemp(find(timesr>t_ini,1):find(timesr>t_f,1));    
    stim(:,jj) = stim(:,jj) - mean(stim(:,jj));
    stim(:,jj) = stim(:,jj)/max(stim(:,jj));
    
    % ORN
    f(:,jj) = msr(find(timesr>t_ini,1):find(timesr>t_f,1));%squeeze(sr(1,1,find(timesr>t_ini,1):find(timesr>t_f,1))); % 
    fsd(:,jj) = sdr(find(timesr>t_ini,1):find(timesr>t_f,1));%zeros(length(f(:,jj)),1); % 
    
    % LN model for ORN response from valve
    fitfun = 1; % type of non-linear static function: 0==(linear fit;) 1==(hill function)
    figon = 1; % 1==(shows plot of LN model and crossprediction)
    reg = 1000; %[10.^[-1:0.5:3]];
    regL2 = 1; % type of regularization: 0==(elastic net); 1==(L2 regularization by me)
    [input_s2n(:,jj), f2_s2n(:,jj), fSTA_s2n(:,jj), fLN_s2n(:,jj), a_s2n(:,jj), D_s2n(:,jj), NR_s2n(:,jj), regopt_s2n(:,jj), lambda_s2n(:,jj)]= build_LNmodel_optreg(time, stim(:,jj), f(:,jj), fsd(:,jj), trials, memory, shift_input, sliding, reg, fitfun, figon, regL2);
   
%     % LN model for ORN response from PID
%     fitfun = 1; % type of non-linear static function: 0==(linear fit;) 1==(hill function)
%     figon = 1; % 1==(shows plot of LN model and crossprediction)
%     reg = 0.005; 
%     regL2 = 1; % type of regularization: 0==(elastic net); 1==(L2 regularization by me)
%     [input_p2n(:,jj), f2_p2n(:,jj), fSTA_p2n(:,jj), fLN_p2n(:,jj), a_p2n(:,jj), D_p2n(:,jj), NR_p2n(:,jj), regopt_p2n(:,jj), lambda_p2n(:,jj)]= build_LNmodel_optreg(time, p(:,jj), f(:,jj), fsd(:,jj), trials, memory, shift_input, sliding, reg, fitfun, figon, regL2);

  
end

%% using p2n_a
delay = 2;

if shift_input>0
    time_temp = time(1:end-shift_input+1);
    p_temp = p(shift_input:end, 2);
    f_temp = f(1:end-shift_input+1, 2);
    fsd_temp = fsd(1:end-shift_input+1, 2);
end
L = length(p_temp);

% with delay
ptest(1:delay) = 0;
ptest(delay+1:L) =p_temp(1:end-delay);

% no delay
% ptest = p_temp;

ptest = ptest*.4 - 0.0037;
ftest = f_temp(memory+1:end);
fsdtest = fsd_temp(memory+1:end);
timetest = time_temp(memory+1:end);
s = zeros(L-memory, memory+1);
for i=1:L-memory
    s(i,:) = ptest(memory+i:-1:i);
end
fSTA = s*a_p2n(1:memory+1,1);
D = D_p2n(:,1);


figure; hold on
plot(fSTA, ftest, '.')
x = 0:0.01:0.1;
y = D(1)*(1./(1+(D(2)./(x-D(4))).^D(3)));
plot(x,y, 'r', 'linewidth', 2)
fLNtest = D(1)*(1./(1+(D(2)./(fSTA-D(4))).^D(3)));
fLNtest(fSTA<D(4))=0; % need this to prevent complex numbers coming out from the hill function
figure; hold on
plot(timetest, ftest, 'k')
plot(timetest, fLNtest, 'r')
Pr = mean(abs(ftest-fLNtest).^2);
Pn = mean(fsdtest.^2);
NR = sqrt(Pr/Pn);

%% using s2n_a
delay = 2;

if shift_input>0
    time_temp = time(1:end-shift_input+1);
    p_temp = stim(shift_input:end, 2);
    f_temp = f(1:end-shift_input+1, 2);
    fsd_temp = fsd(1:end-shift_input+1, 2);
end
L = length(p_temp);

% with delay
ptest(1:delay) = 0;
ptest(delay+1:L) =p_temp(1:end-delay);

% no delay
% ptest = p_temp;

ptest = ptest*.5 + .4; %- 0.0037;
ftest = f_temp(memory+1:end);
fsdtest = fsd_temp(memory+1:end);
timetest = time_temp(memory+1:end);
s = zeros(L-memory, memory+1);
for i=1:L-memory
    s(i,:) = ptest(memory+i:-1:i);
end
fSTA = s*a_s2n(1:memory+1,1);
D = D_s2n(:,1);


figure; hold on
plot(fSTA, ftest, '.')
x = -15:0.5:15;
y = D(1)*(1./(1+(D(2)./(x-D(4))).^D(3)));
plot(x,y, 'r', 'linewidth', 2)
fLNtest = D(1)*(1./(1+(D(2)./(fSTA-D(4))).^D(3)));
fLNtest(fSTA<D(4))=0; % need this to prevent complex numbers coming out from the hill function
figure; hold on
plot(timetest, ftest, 'k')
plot(timetest, fLNtest, 'r')
Pr = mean(abs(ftest-fLNtest).^2);
Pn = mean(fsdtest.^2);
NR = sqrt(Pr/Pn);

%% non-linear part analysis

color = [1 0 0; 0 .75 .75];
figure; 
hold on
for jj=1:num_files
    plot(fSTA_p2n(:,jj), f2_p2n(:,jj), '.', 'color', [.8, .8, .8])
    D = D_p2n(:,jj);
    ffitted = @(x) D(1)*(1./(1+(D(2)./(x-D(4))).^D(3)));
    fplot(ffitted, [min(fSTA_p2n(:,jj)) max(fSTA_p2n(:,jj))]);
end
% plot(fSTA_p2n(:,jj)*.4 - 0.004, f2_p2n(:,jj), '.', 'color', [.8, .8, .8])


%% simulation of filter response to puffs: memory vs time-scale

L = round(4/sliding);
puff = zeros(1,L);
puff(round(2/sliding):round(2.5/sliding)) = 1;
s = zeros(L-memory, memory+1);
for i=1:L-memory
    s(i,:) = puff(memory+i:-1:i);
end

figure; hold on
tau = [0.01:0.03:0.2];
color = jet(length(tau));
for i=1:length(tau);
    a = exp(-(1:memory+1)/(tau(i)/sliding))';
    p_puff = s*a;
%     p_puff = D_s2p(1,i)*p_puff + D_s2p(2,i);
    subplot(2,2,1); hold on
    plot((1:memory+1)*sliding, a, 'color', color(i,:),  'linewidth', 2)
    subplot(2,2,2); hold on
%     plot((1:L)*sliding, puff)
%     subplot(2,2,3); hold on
    plot((memory+1:L)*sliding, p_puff/max(p_puff), 'color', color(i,:),  'linewidth', 2)
end

%% simulation of filter response to puffs: time-scale 

L = round(6/sliding);
tau = 10.^[-0.1:-0.3:-2]/sliding;
amp = max(input_p2n(:,1));
memory2 = memory/2;
a = a_p2n(1:memory2+1,1);
D = D_p2n(:,1);

figure; hold on
subplot(2,3,2)
plot((0:memory2)*sliding, a, 'k', 'linewidth', 2);
xlim([0 .2])
subplot(2,3,3); hold on
ffitted = @(x) D(1)*(1./(1+(D(2)./(x-D(4))).^D(3)));
fplot(ffitted, [D(4) 0.06], 'k');
for j = 1:length(amp)
    color = jet(length(tau));
    for k = 1:length(tau)
        puff = zeros(L,1);
        ind = round(2/sliding):round(2.5/sliding);        
        puff(ind)= amp(j)*(1-exp(-(1:length(ind))/tau(k)));
        ind = round(2.5/sliding):round(6/sliding);
        puff(ind)= max(puff)*exp(-(1:length(ind))/tau(k));
        puff = amp(j)*puff/max(puff);
        s = zeros(L-memory2, memory2+1);
        for i=1:L-memory2
            s(i,:) = puff(memory2+i:-1:i);
        end
        fpuff_lin = s*a;
        subplot(2,3,1); hold on
        plot((1:L)*sliding, puff, 'color', color(k,:), 'linewidth', 2);
        xlim([1.5 5])
        title('stimulus')
        subplot(2,3,6); hold on
        plot((memory2+1:L)*sliding, fpuff_lin, 'color', color(k,:), 'linewidth', 2);
        xlim([1.5 5])
        title('linear response')
       
        fpuff_LN = ffitted(fpuff_lin);
        fpuff_LN(fpuff_lin<D(4))=0;
%         fpuff_LN = fpuff_LN-fpuff_LN(1);
        subplot(2,3,4); hold on
        plot((memory2+1:L)*sliding, fpuff_LN, 'color', color(k,:), 'linewidth', 2);
        xlim([1.5 5])
        subplot(2,3,5); hold on
        plot((memory2+1:L)*sliding, fpuff_LN/max(fpuff_LN), 'color', color(k,:), 'linewidth', 2);
        xlim([1.5 5])
        subplot(2,3,3); hold on
        plot(max(fpuff_lin), ffitted(max(fpuff_lin)), 'o', 'color', color(k,:), 'markerfacecolor', color(k,:))
      
    end
end

%% simulation of filter response to puffs: rescaling

L = round(6/sliding);
tau = 10.^(-2)/sliding;
amp = (0.1:.05:.6)*max(input_p2n(:,1));
memory2 = memory;
a = a_p2n(1:memory2+1,1);
D = D_p2n(:,1);

figure; hold on
subplot(2,3,2)
plot((0:memory2)*sliding, a, 'k', 'linewidth', 2);
xlim([0 .2])
subplot(2,3,3); hold on
ffitted = @(x) D(1)*(1./(1+(D(2)./(x-D(4))).^D(3)));
fplot(ffitted, [0 0.06], 'k');    

color = jet(length(amp));
for j = 1:length(amp)
    for k = 1:length(tau)
        puff = zeros(L,1);
        ind = round(2/sliding):round(2.5/sliding);        
        puff(ind)= amp(j)*(1-exp(-(1:length(ind))/tau(k)));
        ind = round(2.5/sliding):round(6/sliding);
        puff(ind)= max(puff)*exp(-(1:length(ind))/tau(k));
%         puff = max(input_p2n(:,1))*puff/max(puff);
        s = zeros(L-memory2, memory2+1);
        for i=1:L-memory2
            s(i,:) = puff(memory2+i:-1:i);
        end
        fpuff_lin = s*a;
        subplot(2,3,1); hold on
        plot((1:L)*sliding, puff, 'color', color(j,:), 'linewidth', 2);
        xlim([1.5 5])
        title('stimulus')
        subplot(2,3,6); hold on
        plot((memory2+1:L)*sliding, fpuff_lin, 'color', color(j,:), 'linewidth', 2);
        xlim([1.5 5])
        title('linear response')
        
        fpuff_LN = ffitted(fpuff_lin);
        fpuff_LN(fpuff_lin<D(4))=0;
        subplot(2,3,4); hold on
        plot((memory2+1:L)*sliding, fpuff_LN, 'color', color(j,:), 'linewidth', 2);
        xlim([1.5 5])
        subplot(2,3,5); hold on
        plot((memory2+1:L)*sliding,(fpuff_LN-mean(fpuff_LN(1:10)))/max(fpuff_LN-mean(fpuff_LN(1:10))), 'color', color(j,:), 'linewidth', 2);
        xlim([1.5 5])
        subplot(2,3,3); hold on        
        plot(max(fpuff_lin), ffitted(max(fpuff_lin)), 'o', 'color', color(j,:), 'markerfacecolor', color(j,:))
      
    end
end

%% simulation of filter response to puffs: treshold effect

L = round(6/sliding);
tau = 10.^(-2)/sliding;
amp = (0.1:.15:2)*max(input_p2n(:,1));
% tau = 10.^(-.5)/sliding;
% amp = (0.1:.4:5)*max(input_p2n(:,1));

memory2 = memory;
a = a_p2n(1:memory2+1,1);
D = D_p2n(:,1);
D(4) = 0.01;

figure; hold on
subplot(2,3,2)
plot((0:memory2)*sliding, a, 'k', 'linewidth', 2);
xlim([0 .2])
subplot(2,3,3); hold on
ffitted = @(x) D(1)*(1./(1+(D(2)./(x-D(4))).^D(3)));
fplot(ffitted, [0 0.1], 'k');    

color = jet(length(amp));
for j = 1:length(amp)
    for k = 1:length(tau)
        puff = zeros(L,1);
        ind = round(2/sliding):round(2.5/sliding);        
        puff(ind)= amp(j)*(1-exp(-(1:length(ind))/tau(k)));
        ind = round(2.5/sliding):round(6/sliding);
        puff(ind)= max(puff)*exp(-(1:length(ind))/tau(k));
%         puff = max(input_p2n(:,1))*puff/max(puff);
        s = zeros(L-memory2, memory2+1);
        for i=1:L-memory2
            s(i,:) = puff(memory2+i:-1:i);
        end
        fpuff_lin = s*a;
        subplot(2,3,1); hold on
        plot((1:L)*sliding, puff, 'color', color(j,:), 'linewidth', 2);
        xlim([1.5 5])
        title('stimulus')
        subplot(2,3,6); hold on
        plot((memory2+1:L)*sliding, fpuff_lin, 'color', color(j,:), 'linewidth', 2);
        xlim([1.5 5])
        title('linear response')
        
        fpuff_LN = ffitted(fpuff_lin);
        fpuff_LN(fpuff_lin<D(4))=0;
        subplot(2,3,4); hold on
        plot((memory2+1:L)*sliding, fpuff_LN, 'color', color(j,:), 'linewidth', 2);
        xlim([1.5 5])
        subplot(2,3,5); hold on
        plot((memory2+1:L)*sliding,(fpuff_LN-mean(fpuff_LN(1:10)))/max(fpuff_LN-mean(fpuff_LN(1:10))), 'color', color(j,:), 'linewidth', 2);
        xlim([1.5 5])
        subplot(2,3,3); hold on        
        plot(max(fpuff_lin), ffitted(max(fpuff_lin)), 'o', 'color', color(j,:), 'markerfacecolor', color(j,:))
      
    end
end

%% single trials analysis

clear all
files = dir('final_*');
deltat = 10^(-4);

file_list = [1];
jj=0;
for ifile = file_list
    jj = jj+1;
    load(files(ifile).name)
    display(files(ifile).name)
    trialsPID = size(PID,2);
    trials = size(stA,2);
    
    % bin and average
    win = 0.03;
    sliding = 0.003;
    [timesr PIDs sr msr sdr] = bin_traces(PID, stA, deltat, win, sliding);    
    meanPID = mean(squeeze(PIDs(1,:,:)), 1);      
    sdPID = std(squeeze(PIDs(1,:,:)), 1);   
     
%     figure; hold on
%     plot(meanPID, 'k')
%     meanPID = butter_filter(meanPID, .001, 0, sliding, 3, 3);
%     plot(meanPID, 'b')
%     rasterplot(stim_signal, PIDs, meanPID, sdPID/trialsPID,  stA, msr, sdr/trials, deltat, timesr, 10, 20)


    for tr = 1:trials
        t_ini = 9;
        t_f = 35;
        p = squeeze(PIDs(1,tr,find(timesr>t_ini,1):find(timesr>t_f,1)))';
        pid_max(tr,ifile) = max(p);
        p = p - smooth(p, 3/sliding)';
        PIDs(1,tr,find(timesr>t_ini,1):find(timesr>t_f,1)) = p;
        
        % LN model set up
        t_ini = 12;
        t_f = 32;
        memory = round(1/sliding);
        time = timesr(find(timesr>t_ini,1):find(timesr>t_f,1));
        %PID
        p = squeeze(PIDs(1,tr,find(timesr>t_ini,1):find(timesr>t_f,1)));
        pmax = max(p);
        p = p/pmax;
        p = p-mean(p);
        psd = zeros(size(p));
        % ORN
        f = squeeze(sr(1,tr,find(timesr>t_ini,1):find(timesr>t_f,1)))';
        f_max(tr,ifile) = max(smooth(f,10));
        fsd = zeros(size(f));
        % valve
        stim_signal = round(stim_signal);
        [timesr stim] = sliding_average(squeeze(stim_signal(1,1,:)), deltat, win, sliding);
        stim = stim(find(timesr>t_ini,1):find(timesr>t_f,1))';
        stim = stim/max(stim);
        stim = stim - mean(stim);
        
        % LN model for ORN response from valve
%         reg = 1000;
        %     [a_Os, fLN_Os, D_Os, NR_Os, CR_Os] = build_LNmodel_optreg(time, stim, f', fsd', trials, memory, sliding, reg,1, 1);
%         [a(tr,jj,:), fLN(tr,jj,:), D(tr,jj,:), NR(tr,jj,:), CR(tr,jj,:), regopt(tr,jj,:)]= build_LNmodel_optreg(time, stim, f', fsd', trials, memory, sliding, reg, 1, 0);

        %     display(['ORN from valve: NR=' num2str(NR_Os) ' CR='   num2str(CR_Os)]);
        
        %     % LN model for ORN response from PID
            reg = 50;
        %     [a_OP(:,jj), fLN_OP(:,jj), D_OP(:,jj), NR_OP(:,jj), CR_OP(:,jj),regopt(:,jj)] = build_LNmodel_optreg(time, p, f', fsd', trials, memory, sliding, reg, 1, 0);
            [a(:,jj), fLN(:,jj), D(:,jj), NR(:,jj), CR(:,jj), regopt(:,jj)]= build_LNmodel_optreg(time, p, f', fsd', trials, memory, sliding, reg, 1, 0);
        %     display(['ORN from PID: NR=' num2str(NR_OP) ' CR=' num2str(CR_OP)]);
        
        % LN model for PID response from valve
        %     reg = 500;
        %     [a_Ps, fLN_Ps, D_Ps, NR_Ps, CR_Ps] = build_LNmodel_optreg(time, stim, p', psd', trialsPID, memory, sliding, reg, 0, 1);
        %     [a(:,jj), fLN(:,jj), D(:,jj), NR(:,jj), CR(:,jj), regopt(:,jj)] = build_LNmodel_optreg(time, stim, p', psd', trialsPID, memory, sliding, reg, 0, 1);
        %     display(['PID from valve: NR=' num2str(NR_Ps) ' CR=' num2str(CR_Ps)]);
        %}
    end

end

color = [.6 0 0; 0 .5 .5];
figure
for i=1:trials
    subplot(4,5,i); hold on
    for jj = 1:2
        plot((1:memory+1)*sliding, squeeze(a(i,jj,:)), 'linewidth', 2, 'color', color(jj,:))
    end
    xlim([0 .3])
end

figure; hold on
for jj = 1:2
    errorbar_cazzo((1:memory+1)*sliding, mean(squeeze(a(:,jj,:))), std(squeeze(a(:,jj,:)))/sqrt(trials))
end
xlim([0 .3])

color = jet(trials);
figure
for jj = 1:2
    subplot(2,2,jj); hold on
    for i=1:trials
        plot((1:memory+1)*sliding, squeeze(a(i,jj,:)), 'linewidth', 2, 'color', color(i,:))
    end
    xlim([0 .4])
end

delay = [];
for jj = 1:2
    [temp delay(:,jj)] =  max(squeeze(a(:,jj,:))');
end
figure; hold on
plot(delay(:,1)*sliding+ 0.0003*randn(trials,1), delay(:,2)*sliding + 0.0003*randn(trials,1), 'o')
plot([0.06 0.09],[0.06 0.09],'k')
xlim([0.06 0.09])
ylim([0.06 0.09])
[c pvalue] = corr(delay(:,1)*sliding, delay(:,2)*sliding)
    
figure
for jj=1:2
    subplot(2,2,jj); hold on
    plot(pid_max(:,jj),delay(:,jj)*sliding, 'ok')
    [c pvalue] = corr(pid_max(:,jj),delay(:,jj)*sliding)
    fun_fit = fittype('poly1');
    fit_out = fit(pid_max(:,jj), delay(:,jj)*sliding, fun_fit);
    D = coeffvalues(fit_out);
    ffitted = @(x) D(1)*x +D(2);
    fplot(ffitted, [min(pid_max(:,jj)) max(pid_max(:,jj))], 'r');
    subplot(2,2,jj+2); hold on
    plot(f_max(:,jj),delay(:,jj)*sliding, 'ok')
    [c pvalue] = corr(f_max(:,jj),delay(:,jj)*sliding)
    fun_fit = fittype('poly1');
    fit_out = fit(f_max(:,jj), delay(:,jj)*sliding, fun_fit);
    D = coeffvalues(fit_out);
    ffitted = @(x) D(1)*x +D(2);
    fplot(ffitted, [min(f_max(:,jj)) max(f_max(:,jj))], 'r');
end


for jj=1:2
   am(jj,:) = mean(squeeze(a(:,jj,:)));
end

    
    
    
    
    
    
    
    
    
