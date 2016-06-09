%% Testing LFP Models
% In this document, we see if we can write down simple models for the LFP response that show a slowdown in the kinetics and Weber-Fechner gain scalng

pHeader;

% core parameters 
T = 20e3; % total length
pulse_on = 15e3;
pulse_off = 16e3;


%% The Model 
% The first model we consider is as follows: 
% 
% $$ \left(\frac{A+Bs}{1+Cs}\right)\tau\dot{r}=k_{0}+\ln(s)-r $$
% 

%%
% Note that this model can achieve negative responses if s(t) is sufficiently small. This is a little weird, but this can be avoided if we choose $k_{0}$ sufficiently large. 

%%
% First, we check if it does Weber-Fechner gain scaling. We set A = 1, B = 0 and C = 0 to ignore the effects of the timescale for now. 


clear p
p.A = 1;
p.B = 0;
p.C = 0;
p.tau = 1;
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

for i = 1:length(background_levels)
  S = background_levels(i) + zeros(T,1);
  S(pulse_on:pulse_off) = 2*background_levels(i);
  plot(ax(1),S,'Color',c(i,:));

  R = LFPmodel(S,p);
  plot(ax(2),R,'Color',c(i,:))

  temp = R - mean(R(pulse_on-1e3:pulse_on-1));
  plot(ax(3),temp,'Color',c(i,:))

  % compute gain 
  g = max(temp)/(max(S) - min(S));
  plot(ax(4),min(S),g,'+','Color',c(i,:));
end

for i = 1:3
  set(ax(i),'XLim',[length(S) - 10e3 length(S)])
end

set(ax(1),'YScale','log')
ylabel(ax(1),'Stimulus')
ylabel(ax(2),'r(t)')
xlabel(ax(2),'Time')
ylabel(ax(3),'r(t) - r_{background}')
set(ax(3),'YLim',[-.01 .5])
set(ax(4),'XScale','log','YScale','log','XTick',logspace(-2,2,5),'YTick',logspace(-2,2,5))
xlabel(ax(4),'Background S')
ylabel(ax(4),'Gain (\DeltaR / \DeltaS)')


prettyFig('FixLogX','true')

if being_published  
  snapnow  
  delete(gcf)
end

%%
% As expected, this model is capable of changing gain consistent with the Weber-Fechner prediction. Also note that the responses to all fold changes are identical, after removal of background response. 

%%
% Now, we let the timescale effect come into play and see if we can get a slowdown in the kinetics of the response with increasing mean stimulus. We see that as the stimulus becomes very large, the effective timescale scales with B/C, while as the stimulus becomes very small, the effective timescale scales with A. So we set A to be small, and choose B/C to be large. 

clear p
p.A = 1e-3;
p.B = 1;
p.C = 1;
p.tau = 1;
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
end

for i = 1:3
  set(ax(i),'XLim',[length(S) - 10e3 length(S)])
end

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



%% Version Info
%
pFooter;

