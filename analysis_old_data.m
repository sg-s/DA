%% analysis of random stim
% i have rewritten this script, from Carlotta, as a function that accepts a file name
% and spits out the filter, the nonlinear function, the linear fit, the firing rate, 
% and the two predictions
function [model_values] = analysis_old_data(filename)
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

  
    % get the LN model with a non-linear static function
    % fitfun = 1; % type of non-linear static function: 0==(linear fit;) 1==(hill function)
    % figon = 0; % 1==(shows plot of LN model and crossprediction)
    % reg = 50; %[10.^[-1:0.5:3]];
    % regL2 = 1; % type of regularization: 0==(elastic net); 1==(L2 regularization by me)
    % [model_values(1)]= build_LNmodel_optreg(time, p(:,jj), f(:,jj), fsd(:,jj), trials, memory, shift_input, sliding, reg, fitfun, figon, regL2);

    % get the LN model with a linear function at the end
    fitfun = 0; % type of non-linear static function: 0==(linear fit;) 1==(hill function)
    figon = 0; % 1==(shows plot of LN model and crossprediction)
    reg = 50; %[10.^[-1:0.5:3]];
    regL2 = 1; % type of regularization: 0==(elastic net); 1==(L2 regularization by me)

    [model_values(2)]= build_LNmodel_optreg(time, p(:,jj), f(:,jj), fsd(:,jj), trials, memory, shift_input, sliding, reg, fitfun, figon, regL2);

    % add some more useful data
    % add the PID and stimulus

    model_values(1).PID =p(find(time==min(model_values(1).time)):find(time==max(model_values(1).time)));
    % carlotta's code mangles the stimulus in some weird way. so need to raclcualte
    stim = (squeeze(stim_signal));
    if ~isvector(stim)
        stim = mean(stim);
    end   
    stim=stim(min(model_values(1).time)*1e4:max(model_values(1).time)*1e4);
    model_values(1).stim = stim(1:floor(sliding/1e-4):end);

    % calcualte on and off durations
    stim(stim>0.5) = 1;
    stim(stim<=0.5) = 0;
    d = diff(stim);

    ons = find(d>0.5);
    offs = find(d<-0.5);
    if offs(1) < ons(1)
        offs(1) = [];
    end
    if length(ons) > length(offs) && ons(end) > offs(end)
        ons(end) = [];
    end
    if length(ons) ~= length(offs)
        beep
        keyboard
    end
    model_values(1).on_durations = (offs-ons)*1e-4;
    ons(1) = [];
    offs(end) = [];
    model_values(1).off_durations = (ons-offs)*1e-4;

  
end


return
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

    
    
    
    
    
    
