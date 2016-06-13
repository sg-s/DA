
%% The Nagel-Wilson Model
% In this document we try to make sense of the Nagel Wilson Model and see if it explains gain scaling and LFP slowdown.

pHeader;

% core parameters 
T = 20e3; % total length
pulse_on = 15e3;
pulse_off = 16e3;

%% 1: Using the Nagel-Wilson Model as is
% In this section, we use the Nagel-Wilson model as described in the paper. This has 6 ODEs, and the output of the model is the channel open probability. In the following figure, I plot the response of this model to a test pulse of stimulus. 

clear p
p.theta = 100;
p.   ka = 0.001;
p.   kb = 1;
p.   sa = 10;
p.   sb = 100;
p.   k0 = 500;
p.    A = 0.8000;
p.    B = 0.6000;
p.   ko = 500;
p.   kc = 10;

S = .1+zeros(T,1);
S(pulse_on:pulse_off) = 1;
[C,D,R] = NagelWilsonIntegrate(S,p);

figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
subplot(2,2,1), hold on
plot(S)
ylabel('Stimulus')
subplot(2,2,2), hold on
plot(R)
legend({'R','R*','OR','OR*'},'Location','northwest')
ylabel('Receptor Species')
subplot(2,2,3), hold on
plot(C)
ylabel('C_{open}')
subplot(2,2,4), hold on
plot(D)
ylabel('Diffusible factor')

prettyFig();

if being_published
  snapnow
  delete(gcf)
end

%%
% We see that there is a fatal flaw in the model: there is no conservation equation for the four reaction species. Thus, they are allowed to grow without bound, which is obviously wrong. 

%% 2. Nagel Wilson Model with conservation equation
% In this section, we add a conservation equation to the model:
% 
% $$ OR=R_{total}-R-R*-OR*$$
%

%%
% Doing so eliminates one of the ODEs. The ODEs now are:
%
% $$ \dot{R}=-R(k_{a}s_{a}+Ok_{b}s_{b})+R*s_{a}+(R_{total}-R-R*-OR*)s_{b} $$
%

%%
% and so on.

%%
% In the following figure, we plot the responses to a pulse of odorant.

S = .1+zeros(T,1);
S(pulse_on:pulse_off) = .2;
[C,D,R] = NagelWilsonIntegrate2(S,p);

figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
subplot(2,2,1), hold on
plot(S)
set(gca,'XLim',[length(S)-10e3 length(S)])
ylabel('Stimulus')
subplot(2,2,2), hold on
plot(R)
set(gca,'XLim',[length(S)-10e3 length(S)])
legend({'R','R*','OR','OR*'})
ylabel('Receptor Species')
subplot(2,2,3), hold on
plot(C)
set(gca,'XLim',[length(S)-10e3 length(S)])
ylabel('C_{open}')
subplot(2,2,4), hold on
plot(D)
set(gca,'XLim',[length(S)-10e3 length(S)])
ylabel('Diffusible factor')

prettyFig();

if being_published
  snapnow
  delete(gcf)
end

%%
% OK, this is somewhat reasonable. Now, we check to see if this model predicts a slowdown in response with increased stimulus. We do this by stimulating the model with pulses on top of increasing backgrounds. We also check to see if the responses follow Weber's Law.

clear ax
figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
for i = 1:6
  ax(i) = subplot(2,3,i); hold on
end

background_levels = logspace(-2,2,10);
c = parula(length(background_levels)+1);

% plot 1/S guide lines
guide_A = logspace(-4,-1,4);
for i = 1:length(guide_A)
  plot(ax(5),background_levels,guide_A(i)./background_levels,'r--') 
  plot(ax(6),background_levels,guide_A(i)./background_levels,'r--')
end

for i = 1:length(background_levels)
  S = background_levels(i) + zeros(T,1);
  S(pulse_on:pulse_off) = 2*background_levels(i);
  plot(ax(1),S,'Color',c(i,:));

  [C,D,R] = NagelWilsonIntegrate2(S,p);
  plot(ax(2),R(:,2)+R(:,4),'Color',c(i,:))

  plot(ax(3),C,'Color',c(i,:))

  plot(ax(4),D,'Color',c(i,:))

  % compute gain at receptors 
  y = R(:,2) + R(:,4);
  delta_r = max(y(pulse_on:pulse_off)) - y(pulse_on-1);
  g = delta_r/(max(S) - min(S));
  plot(ax(5),min(S),g,'+','Color',c(i,:),'MarkerSize',24);

  % plot total gain
  delta_r = max(C(pulse_on:pulse_off)) - C(pulse_on-1);
  g = delta_r/(max(S) - min(S));
  plot(ax(6),min(S),g,'+','Color',c(i,:),'MarkerSize',24);
end

for i = 1:4
  set(ax(i),'XLim',[length(S) - 10e3 length(S)])
end

set(ax(1),'YScale','log')
ylabel(ax(1),'Stimulus')
ylabel(ax(2),'R* + OR*')
ylabel(ax(3),'C_{open}')
ylabel(ax(4),'Diffusible factor')

title(ax(5),'Receptor Gain')
ylabel(ax(5),'\DeltaR/\DeltaS')
set(ax(5),'YScale','log','XScale','log')
xlabel(ax(5),'Background Stim')

title(ax(6),'Channel Gain')
ylabel(ax(6),'\DeltaR/\DeltaS')
set(ax(6),'YScale','log','XScale','log')
xlabel(ax(6),'Background Stim')

prettyFig();

if being_published
  snapnow
  delete(gcf)
end

%%
% This model doesn't follow Weber's Law, at least for these parameters. Gain appears to fall off faster than expected. 


%%
% That looks interesting. What if we rescale the channel open probabilities to look more closely at the kinetics? 


clear ax
figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
for i = 1:4
  ax(i) = subplot(2,2,i); hold on
end

for i = 1:length(background_levels)
  S = background_levels(i) + zeros(T,1);
  S(pulse_on:pulse_off) = 2*background_levels(i);
  plot(ax(1),S,'Color',c(i,:));

  [C,D,R] = NagelWilsonIntegrate2(S,p);

  R = R(:,2) + R(:,4);
  R = R - mean(R(pulse_on-1e3:pulse_on));
  R = R/max(R(pulse_on:pulse_off));
  plot(ax(2),R,'Color',c(i,:))

  C = C - mean(C(pulse_on-1e3:pulse_on));
  C = C/max(C(pulse_on:pulse_off));

  plot(ax(3),C,'Color',c(i,:));

  D = D - mean(D(pulse_on-1e3:pulse_on));
  D = D/max(D(pulse_on:pulse_off));

  plot(ax(4),D,'Color',c(i,:))
end

ylabel(ax(1),'Stimulus')
ylabel(ax(2),'R* + OR* (rescaled)')
ylabel(ax(3),'C_{open} (rescaled)')
ylabel(ax(4),'Diffusible Factor (rescaled)')
set(ax(1),'XLim',[length(S) - 10e3 length(S)],'YScale','log')
set(ax(2),'XLim',[length(S) - 10e3 length(S)],'YLim',[-.1 1.1])
set(ax(3),'XLim',[length(S) - 10e3 length(S)],'YLim',[-.5 1.1])
set(ax(4),'XLim',[length(S) - 10e3 length(S)],'YLim',[-.1 1.1])

prettyFig();

if being_published
  snapnow
  delete(gcf)
end


%%
% The time course of the response seems to change, but it doesn't look like the kinetics seem to slow down with increasing background. If anything, the time to peak decreases with increasing background stimulus. 





%% Version Info
%
pFooter;
