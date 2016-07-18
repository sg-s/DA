

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
p.    B = 20.7656;

%% Fitting this model to simplified DA model simulations 
% In the first section, we create some synthetic data by running the reduced DA model (without the ODE) on some Gaussian white noise, and then adding some noise. We then attempt to recover the filters of the DA model. The following figure shows the DA model filters, together with the recovered filters. First, we do this without whitening the stimulus (as we don't need to).

%%
% We can directly compute the filters from:
%
% $$ K=R/\hat{S} $$
%

%%
% where $K$ is a $1x2N$ vector, containing both filters one after the other, and $R$ is a $1xT$ vector of the responses, and $\hat{S}$ is a $2NxT$ matrix formed from concatenating the time-shifted stimulus and the time-shifted stimulus times the response. 

RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1984)); 
S = randn(20e3,1);
[R,~,~,Ky,Kz] = DAModelv2(S,p);
R = R + .1*randn(length(R),1);

[K0,K1] = fitFirstOrderModel(S(5e3:end),R(5e3:end),'filter_length',1e3,'whiten',false,'reg',0);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
plot(K0,'r')
plot(Ky*p.A,'k')

xlabel('Filter lag (a.u.)')
legend({'K_0','\alpha K_y'});
subplot(1,2,2); hold on
plot(-K1,'r')
plot(Kz*p.B,'k'); hold on
xlabel('Filter lag (a.u.)')
legend({'-K_1','\beta K_z'});

suptitle('Gaussian white inputs, +noise, No whitening')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% Now we repeat the same thing, but whiten the stimulus and back out filters. These filters are given by:
%
% $$ K=C\setminus(\hat{S} *R) $$
%

%%
% where $C$ is the covariance matrix of the time-shifted stimulus. 


[K0,K1] = fitFirstOrderModel(S(5e3:end),R(5e3:end),'filter_length',1e3,'whiten',true,'reg',0);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
plot(K0,'r')
plot(Ky*p.A,'k')

xlabel('Filter lag (a.u.)')
legend({'K_0','\alpha K_y'});
subplot(1,2,2); hold on
plot(-K1,'r')
plot(Kz*p.B,'k'); hold on
xlabel('Filter lag (a.u.)')
legend({'-K_1','\beta K_z'});

suptitle('Gaussian inputs, whitening')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% So that works too. Now we do the same thing, but the stimulus is now correlated over a short timescale. This time, we attempt to regularise the solution by diagonalising the covarince matrix over each filter block. 

tc = 50;
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',3)); 
S = filtfilt(ones(tc,1),tc,randn(20e3,1));
[R,~,~,Ky,Kz] = DAModelv2(S,p);
R = R + .1*randn(length(R),1);

[K0,K1] = fitFirstOrderModel(S(5e3:end),R(5e3:end),'filter_length',1e3,'whiten',true,'reg',1e-4);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
plot(K0,'r')
plot(Ky*p.A,'k')

xlabel('Filter lag (a.u.)')
legend({'K_0','\alpha K_y'});
subplot(1,2,2); hold on
plot(-K1,'r')
plot(Kz*p.B,'k'); hold on
xlabel('Filter lag (a.u.)')
legend({'-K_1','\beta K_z'});

suptitle('correlated inputs, small regularisation')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%% Version Info
%
pFooter;


