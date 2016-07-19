

pHeader;

%% nth-Order Models
% In this document, we write down a class of models where successive orders of the response depend on various linear filters convolved with the stimulus. The model reduces to the linear model for zeroth order, the simplified DA model for the 1st order, and the full DA model for the 2nd order. This is what the equations look like
%
%
% $$ 0=K_{0}\otimes s(t)+r(t) $$
%
% $$ 0=K_{0}\otimes s(t)+(\alpha_{1}+K_{1}\otimes s(t))\cdot r(t) $$
%
% $$ 0=K_{0}\otimes s(t)+(\alpha_{1}+K_{1}\otimes s(t))\cdot r(t)+(\alpha_{2}+K_{2}\otimes s(t))\cdot\frac{dr(t)}{dt} $$ 
%


clear p
p.   s0 = 0;
p.  n_z = 2;
p.tau_z = 150;
p.  n_y = 2;
p.tau_y = 30;
p.    C = 0.5000;
p.    A = 170.4375;
p.    B = 2.7656;

%% Fitting this model to simplified DA model simulations 
% In the first section, we create some synthetic data by running the reduced DA model (without the ODE) on some Gaussian white noise, and then adding some noise. We then attempt to recover the filters of the DA model. The following figure shows the DA model filters, together with the recovered filters. First, we do this without whitening the stimulus (as we don't need to).

%%
% We can directly compute the filters from:
%
% $$ K=R/\hat{S} $$
%

%%
% It is more efficient to compute the filters using:
%
% $$ K=C\setminus(\hat{S} *R) $$
%

%%
% where C is the covariance matrix of the shifted stimulus. 

%%
% where $K$ is a $1\times2N$ vector, containing both filters one after the other, and $R$ is a $1 \times T$ vector of the responses, and $\hat{S}$ is a $2N \times T$ matrix formed from concatenating the time-shifted stimulus and the time-shifted stimulus times the response. 

RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1984)); 
S = randn(20e3,1);
[R,~,~,Ky,Kz] = DAModelv2(S,p);
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1984));
R = R + .5*randn(length(R),1);

filterset = fitNOrderModel(S(5e3:end),R(5e3:end),'filter_length',1e3,'reg',0,'offset',100);
R_pred = firstOrderModel(S,filterset);


figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on
subplot(1,3,1); hold on
plot(filterset.time,filterset.K0,'r')
plot(Ky*p.A,'k')

xlabel('Filter lag (a.u.)')
legend({'K_0','\alpha K_y'});
subplot(1,3,2); hold on
plot(filterset.time,filterset.K1,'r')
plot(Kz*p.B,'k'); hold on
xlabel('Filter lag (a.u.)')
legend({'K_1','\beta K_z'});

subplot(1,3,3); hold on
plot(R_pred(1e3:end),R(1e3:end),'.')
legend(['r^2 = ' oval(rsquare(R_pred,R))],'Location','southeast')
xlabel('Prediction')
ylabel('Data')

suptitle('Gaussian white inputs, +noise, No regularisation')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end



%%
% Now I do the same thing, but the stimulus is now correlated over a short timescale. Regularization is attempted by adding white noise to the stimulus and then extracting filters. The regularisation factor then is the scale of the white noise added to the stimulus. The best regularisation factor is found through cross-validation. In addition, I had to employ a number of hacks to get this procedure working so that decent looking filters that predict the response well were obtained. They are:
% 
% # I had to make sure that the extracted filters were long enough. If they're not, the procedure breaks down (see the next steps)
% # I allowed for an offset. The reason is because I had to trim the filter to kill some ringing at the ends (see next)
% # I had to trim the ends to kill some ringing at the ends. 
% # I smoothed the filter with some small timescale to prevent filters from being dominated by high-frequency components. I found that smoothing, then trimming gave better results than trimming and then smoothing (which is weird)
% 

%%
% Despite these hacks, the procedure works well in that no tuning is required to get the filters out: it runs automatically on the data, and generates very good approximations of the filters. See below. 

tc = 50;
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1984));
S = filtfilt(ones(tc,1),tc,randn(20e3,1));
[R,~,~,Ky,Kz] = DAModelv2(S,p);
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1984));
R = R + .1*randn(length(R),1);

[filterset, diagnostics] = fitNOrderModel(S(5e3:end),R(5e3:end),'filter_length',1200,'reg','best','filter_low_pass',50);

figure('outerposition',[0 0 1500 501],'PaperUnits','points','PaperSize',[1500 501]); hold on
subplot(1,3,1); hold on
plot(filterset.time,filterset.K0,'r')
plot(Ky*p.A,'k')
xlabel('Filter lag (a.u.)')
legend({'K_0','\alpha K_y'});

subplot(1,3,2); hold on
plot(filterset.time,filterset.K1,'r')
plot(Kz*p.B,'k'); hold on
xlabel('Filter lag (a.u.)')
legend({'K_1','\beta K_z'});

subplot(1,3,3); hold on
plot(diagnostics.reg_vec,diagnostics.r2,'k+')
set(gca,'XScale','log','YScale','log','YLim',[.1 1],'XTick',[])
xlabel('Regularisation factor')
ylabel('r^2')

suptitle('correlated inputs, adding white noise to inputs')

prettyFig('FixLogX',true)

if being_published	
	snapnow	
	delete(gcf)
end

%%
% Note that the filters "bleed" a little before 0 -- which is a consequence of the smoothing that I have to do. 

%% Real Data
% Now, I attempt to run it on some real data. First, I run it on DA model responses to naturalistic stimuli. Note that I fit this DA model to the ORN response to this stimulus, and it does a very good job at it. 


load(getPath(dataManager,'5c7dacc5b42ff0eebb980d80fec120c3'),'data','spikes')
PID = data(2).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(2).A;

% A spikes --> firing rate
fA = spiketimes2f(all_spikes,time);

tA = 1e-3*(1:length(fA));
PID2 = fA;
for i = 1:width(PID2)
	PID2(:,i) = interp1(time,PID(i,:),tA);
end
PID = PID2; clear PID2
% some minor cleaning up
PID(end,:) = PID(end-1,:); 

% remove the baseline from the PID, and remember the error
PID = PID - mean(mean(PID(1:5e3,:)));
S = nanmean(PID,2);

[R,~,~,Ky,Kz] = DAModelv2(S,p);
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1984)); 
R = R + .1*randn(length(R),1);

filterset = fitNOrderModel(S(5e3:end),R(5e3:end),'filter_length',1200,'reg','best','filter_low_pass',50);

figure('outerposition',[0 0 1500 801],'PaperUnits','points','PaperSize',[1500 801]); hold on
subplot(2,3,1:3); hold on
plot(R,'k') 
R_pred = firstOrderModel(S,filterset);
plot(R_pred,'r')
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
legend({'DA Model Response','1st order model'})
text(1e4,100,['r^2=' oval(rsquare(R,R_pred))])
set(gca,'XLim',[0 70e3])

subplot(2,3,4); hold on
plot(filterset.time,filterset.K0,'r')
plot(Ky*p.A,'k')
xlabel('Filter lag (a.u.)')
legend({'K_0','\alpha K_y'});

subplot(2,3,5); hold on
plot(filterset.time,filterset.K1,'r')
plot(Kz*p.B,'k'); hold on
xlabel('Filter lag (a.u.)')
legend({'K_1','\beta K_z'});

subplot(2,3,6); hold on
plot(R_pred(1:10:end),R(1:10:end),'k.')
xlabel('Prediction')
ylabel('DA Model response')

suptitle('naturalistic inputs, DA model responses')

prettyFig('FixLogX',true)

if being_published	
	snapnow	
	delete(gcf)
end

%%
% OK, that looks good. Now we fit the model to the real data (ORN responses)

R = mean(fA,2);

filterset = fitNOrderModel(S(5e3:end),R(5e3:end),'filter_length',1200,'reg','best','filter_low_pass',20,'cross_validate',false);

figure('outerposition',[0 0 1500 801],'PaperUnits','points','PaperSize',[1500 801]); hold on
subplot(2,3,1:3); hold on
plot(R,'k') 
R_pred = firstOrderModel(S,filterset);
plot(R_pred,'r')
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
legend({'ab3A','1st order model'})
text(1e4,100,['r^2=' oval(rsquare(R,R_pred))])
set(gca,'XLim',[0 70e3])

subplot(2,3,4); hold on
plot(filterset.time,filterset.K0,'r')
xlabel('Filter lag (a.u.)')
legend({'K_0'});

subplot(2,3,5); hold on
plot(filterset.time,filterset.K1,'r')
xlabel('Filter lag (a.u.)')
legend({'K_1'});

subplot(2,3,6); hold on
plot(R_pred(1:10:end),R(1:10:end),'k.')
xlabel('Prediction')
ylabel('ab3A firing rate (Hz)')

suptitle('naturalistic inputs + real ORN data')

prettyFig('FixLogX',true)

if being_published	
	snapnow	
	delete(gcf)
end

%%
% Now we fit it to another naturalistic stimulus dataset. 

load(getPath(dataManager,'aeb361c027b71938021c12a6a12a85cd'),'-mat');
example_orn = 4;
S = nanmean(od(example_orn).stimulus,2);
R = nanmean(od(example_orn).firing_rate,2);
S = S - nanmean(S(1:5e3));

filterset = fitNOrderModel(S(5e3:end-5e3),R(5e3:end-5e3),'filter_length',700,'reg',.11,'filter_low_pass',10,'cross_validate',false,'offset',150,'left_trim',50,'right_trim',50);

figure('outerposition',[0 0 1500 801],'PaperUnits','points','PaperSize',[1500 801]); hold on
subplot(2,3,1:3); hold on
plot(R,'k') 
R_pred = firstOrderModel(S,filterset);
plot(R_pred,'r')
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
legend({'ab3A','1st order model'})
text(3e3,100,['r^2=' oval(rsquare(R,R_pred))])
set(gca,'XLim',[0 70e3])

subplot(2,3,4); hold on
plot(filterset.time,filterset.K0,'r')
xlabel('Filter lag (a.u.)')
legend({'K_0'});

subplot(2,3,5); hold on
plot(filterset.time,filterset.K1,'r')
xlabel('Filter lag (a.u.)')
legend({'K_1'});

subplot(2,3,6); hold on
plot(R_pred(1:10:end),R(1:10:end),'k.')
xlabel('Prediction')
ylabel('ab3A firing rate (Hz)')

suptitle('naturalistic inputs + real ORN data')

prettyFig('FixLogX',true)

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


