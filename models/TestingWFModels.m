%% Testing Models
% In this document, we see if we can write down simple models for the LFP response that show a slowdown in the kinetics and Weber-Fechner gain scaling

pHeader;

% core parameters 
T = 20e3; % total length
pulse_on = 15e3;
pulse_off = 16e3;


%% 1. LFP Model
% The first model we consider is as follows: 
% 
% $$ \left(\frac{A+Bs}{1+Cs}\right)\tau\dot{r}=k_{0}+\ln(s)-r $$
% 

%%
% Note that this model can achieve negative responses if s(t) is sufficiently small. This is a little weird, but this can be avoided if we choose $k_{0}$ sufficiently large. 


%%
% Now, we let the timescale effect come into play and see if we can get a slowdown in the kinetics of the response with increasing mean stimulus. We see that as the stimulus becomes very large, the effective timescale scales with B/C, while as the stimulus becomes very small, the effective timescale scales with A. So we set A to be small, and choose B/C to be large. 

clear p
p.A = 1e3;
p.B = 1;
p.tau = 1e-3;
p.ko = 1;


clear ax
figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
for i = 1:4
  ax(i) = subplot(2,2,i); hold on
end

background_levels = logspace(-2,2,10);
c = parula(length(background_levels)+1);

% make sure we don't get negative responses
p.ko = -log(background_levels(1))*1.1;

all_gain = NaN*background_levels;
all_mu = NaN*background_levels;
for i = 1:length(background_levels)
  S = background_levels(i) + zeros(T,1);
  S(pulse_on:pulse_off) = 2*background_levels(i);
  plot(ax(1),S,'Color',c(i,:));

  R = LFPmodel(S,p);
  plot(ax(2),R,'Color',c(i,:))

  temp = R - mean(R(pulse_on-1e3:pulse_on-1));

  plot(ax(3),temp/max(temp),'Color',c(i,:))


  % compute gain 
  g = max(temp)/(max(S) - min(S));
  plot(ax(4),min(S),g,'+','Color',c(i,:));
  all_gain(i) = g;
  all_mu(i) = max(S) - min(S);
end

for i = 1:3
  set(ax(i),'XLim',[length(S) - 10e3 length(S)])
end

% draw a weber's fit to this
x = all_mu;
y = all_gain;
options = fitoptions(fittype('poly1'));
options.Lower = [-1 -Inf];
options.Upper = [-1 Inf];
cf_temp = fit(log(x(:)),log(y(:)),'poly1',options);
cf = fit(x(:),y(:),'power1');
warning off
cf.a = exp(cf_temp.p2); cf.b = -1;
warning on
plot(ax(4),sort(x),cf(sort(x)),'r')

set(ax(1),'YScale','log')
ylabel(ax(1),'Stimulus')
ylabel(ax(2),'r(t)')
xlabel(ax(2),'Time')
ylabel(ax(3),'r(t) (rescaled)')
set(ax(3),'YLim',[-.1 1])
set(ax(4),'XScale','log','YScale','log','XTick',logspace(-2,2,5),'YTick',logspace(-2,2,5))
xlabel(ax(4),'Background S')
ylabel(ax(4),'Gain (\DeltaR / \DeltaS)')


prettyFig('FixLogX','true')

if being_published  
  snapnow  
  delete(gcf)
end

%% 
% It looks like the kinetics are slowing down, as expected. 

%% 2. Schulze-Louis Model 
% Here we test models from Mathieu Louis' eLife paper. I tried using the parameters in their paper, but those parameters didn't yield Weber-Fechner scaling. Fitting their model to our data was very hard, as it's very expensive to solve the model. Instead, I played with the parameters a bit, and found a set that gave Weber-Fechner gain scaling: 


clear p
p.   a1 = 2;
p.   a2 = 1;
p.   a3 = 0;
p.   b1 = 100;
p.   b2 = 8.6300;
p.   b3 = 2.4000;
p.   b4 = 0.5000;
p.   b5 = 0.5000;
p.theta = 1;
p.    n = 2;

clear ax
figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
for i = 1:4
  ax(i) = subplot(2,2,i); hold on
end

background_levels = logspace(0,2,10);
c = parula(length(background_levels)+1);

all_gain = NaN*background_levels;
all_mu = NaN*background_levels;
for i = 1:length(background_levels)
  S = background_levels(i) + zeros(T,1);
  S(pulse_on:pulse_off) = 2*background_levels(i);
  plot(ax(1),S,'Color',c(i,:));

  R = SchulzeLouisModel(S,p);
  plot(ax(2),R,'Color',c(i,:))

  temp = R - mean(R(pulse_on-1e3:pulse_on-1));

  plot(ax(3),temp/max(temp(pulse_on:pulse_off)),'Color',c(i,:))

  % compute gain 
  g = max(temp(pulse_on:pulse_off))/(max(S) - min(S));
  plot(ax(4),min(S),g,'+','Color',c(i,:));
  all_gain(i) = g;
  all_mu(i) = max(S) - min(S);
end

for i = 1:3
  set(ax(i),'XLim',[length(S) - 10e3 length(S)])
end

% draw a weber's fit to this
x = all_mu;
y = all_gain;
options = fitoptions(fittype('poly1'));
options.Lower = [-1 -Inf];
options.Upper = [-1 Inf];
cf_temp = fit(log(x(:)),log(y(:)),'poly1',options);
cf = fit(x(:),y(:),'power1');
warning off
cf.a = exp(cf_temp.p2); cf.b = -1;
warning on
plot(ax(4),sort(x),cf(sort(x)),'r')


set(ax(1),'YScale','log')
ylabel(ax(1),'Stimulus')
ylabel(ax(2),'r(t)')
xlabel(ax(2),'Time')
ylabel(ax(3),'r(t) (rescaled)')
set(ax(3),'YLim',[-.1 1])
set(ax(4),'XScale','log','YScale','log','XTick',logspace(-2,2,5),'YTick',logspace(-2,2,5))
xlabel(ax(4),'Background S')
ylabel(ax(4),'Gain (\DeltaR / \DeltaS)')


prettyFig('FixLogX','true')

if being_published  
  snapnow  
  delete(gcf)
end




%% Version Info
%
pFooter;

