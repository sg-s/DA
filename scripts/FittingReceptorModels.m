% FittingReceptorModels
% 
% created by Srinivas Gorur-Shandilya at 3:31 , 26 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


%% Fitting Simple Receptor Models to Data
% In this document, we look at very simple receptor-based and phenomenological models and see if they can account for 
% 
% # gain changes with mean stimulus 
% # response speedups with increasing mean stimulus. 
% # gain changes with stimulus contrast 
% # responses to triangle waveforms 
% 

% add homebrew path
path1 = getenv('PATH');
if isempty(strfind(path1,':/usr/local/bin'))
    path1 = [path1 ':/usr/local/bin'];
end
setenv('PATH', path1);

% this code determines if this function is being called
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
		unix(['tag -a publish-failed ',which(mfilename)]);
		unix(['tag -r published ',which(mfilename)]);
	end
end
tic


load('/code/da/data/MSG_per_neuron')

% convert data into a nicer format
clear data
data = struct;
data(1).stimulus = mean2([MSG_data(1,:).stim]);
data(1).response = mean2([MSG_data(1,:).resp]);
for i = 2:8
	%data(i).stimulus = mean2([MSG_data(i,:).stim]);
	%data(i).response = mean2([MSG_data(i,:).resp]);
	%data(i).response(1:5e3) = NaN;
	%data(i).stimulus = mean(mean([MSG_data(i,:).stim])) + mean2([MSG_data(1,:).stim]);
	data(i).stimulus = i*.1 + mean2([MSG_data(1,:).stim]);
end

%% Minimal Receptor Model
% In this model, we consider a receptor that can bind and unbind with the odorant signal. The bound fraction is converted to the firing rate using a Hill function. That's it. Here is the ODE governing this sytem:
%
% $$\dot{b} = \alpha (1-b) s(t) - \alpha \theta_{b}b $$
% $$ r = A\frac{b^{2}}{b^2 + K_d^2} $$


%%
% The reason we don't simply add on a LN model to the receptor binding is as follows:

%%
% If we considered our model to be simply receptor binding, and regressed a LN model to the model prediction and the data, there are two problems. First, we would always converge to a globally attractive solution where receptor binding rates would vanish, and we would simply get back a LN model. Second, we can't use this model to generate responses to stimuli, as the LN model is not explicitly accounted for.

%%
% If, on the other hand, we parameterized the LN model, and fit the whole thing to the data, we would never converge to a solution. This is because the receptor ODE and the filter essentially trade-off one for the other, meaning that there infinitely many solutions between various contributions of the two. One way of thinking of this is that we are fitting a model with two parameters where the two parameters only ever occur in the model as a product. So there is no solution where the individual values of the two parameters matter -- only the product matters. 

%%
% Anyway, this simple receptor model can account for response speedups with increasing stimulus (as the response rate scales with the stimulus). But the form of how the gain varies with the stimulus is all wrong. 

clear p
p.    r_b = 0.0036;
p.theta_b = 2.9258;
p. hill_A = 561.5096;
p. hill_K = 0.9528;

characteriseModel(@simpleReceptorModel,p,data);

prettyFig('fs=14;')

if being_published
	snapnow
	delete(gcf)
end


%% Receptors with -ve Feedback
% We now add a diffusible factor that increases the rate of receptor unbinding. We do this because this is the only consistent way of getting a negative feedback with increasing speeds on increasing stimulus. If instead we decreased the binding rate (e.g., by dividing by the diffusible factor), then responses would get slower with increasing stimulus. 

%%
% The ODE governing the diffusible factor is the simplest possible: diffusible factor grows with the stimulus, and decays linearly. We see that this extremely simple model (2 linear ODEs, 6 parameters) can show a non-trivial gain scaling that varies as a power law (close to -1), and that shows response speedups with increasing stimulus. 
%
% $$ \dot{b} = \alpha (1-b) s(t) - \alpha \theta_{b}b d $$
% $$ \dot{d} = \beta s(t) - \beta \theta_{d} d  $$
% $$ r = A\frac{b^{2}}{b^2 + K_d^2} $$ 

clear p
p.    r_b = 0.0013;
p.    r_d = 0.0134;
p.theta_b = 1.0000;
p.theta_d = 0.0820;
p. hill_A = 1.1878e+03;
p. hill_K = 0.5805;

characteriseModel(@simpleReceptorModelv2,p,data);

prettyFig('fs=14;')

if being_published
	snapnow
	delete(gcf)
end

%% 
% For completeness, we consider a model where the binding rate decreases with increasing odor stimulus. As we predicted, this model no longer shows response speedups with increasing stimuli. In fact, because the diffusible factor scales with the stimulus, the response timescale is stimulus-independent. 
%
% $$ \dot{b} =\frac{\alpha (1-b) s(t) }{1+d} - \alpha \theta_{b}b $$
% $$ \dot{d} = \beta s(t) - \beta \theta_{d} d  $$
% $$ r = A\frac{b^{2}}{b^2 + K_d^2} $$ 

clear p
p.    r_b = 0.0099;
p.    r_d = 0.0156;
p.theta_b = 1.0239;
p.theta_d = 0.0820;
p. hill_A = 1.5150e+03;
p. hill_K = 0.5772;

characteriseModel(@simpleReceptorModelv3,p,data);

prettyFig('fs=14;')

if being_published
	snapnow
	delete(gcf)
end

%% Weber-Derived Model
% In this section, we explicitly force the diffusible factor to depend only on the binding fraction, instead of the stimulus, as the previous assumption was a bit unrealistic. Instead, we let the rate of diffusible factor increase depend generally on the bound fraction, and impose the Weber-Fechner Law, and derive the relationship. 
% 
% $$ \dot{b} =\alpha (1-b) s(t) - \alpha \theta_{b}b d $$
% $$ \dot{d} = \beta \frac{1-b}{b}e^{b/k} - \beta \theta_{d} d $$
% $$ r = A*b + C $$ 

%%
% This model has a problem when $b=0$, and the diffusible factor explodes (you can see this happen in some cases).

clear p
p.    r_b =  0.0204;
p.    r_d =  0.0239;
p.theta_b =  4.5938;
p.theta_d =  13.3594;
p.      k =  0.6240;
p.      A =  59.7812;
p.     R0 =  -3.4307;
characteriseModel(@longReceptorModel,p,data);

prettyFig('fs=14;')


if being_published
	snapnow
	delete(gcf)
end


%% DA Model
% In this section, we fit the DA model to the data and see what it can do. The DA Model is described by:
% $$ r(t) = \frac{\alpha K_{y}\otimes s(t)}{1+\beta K_{z}\otimes{s(t)}} $$


clear p
p.tau_z =  128.6250;
p.tau_y =  27.3313;
p.  n_y =  2;
p.  n_z =  2;
p.    A =  155.0387;
p.    B =  4.8162;
p.    C =  0.0454;
p.   s0 =  0;
characteriseModel(@DAModelv2,p,data);

prettyFig('fs=14;')


if being_published
	snapnow
	delete(gcf)
end

%% DA Model with contrast adaptation
% In this section, we modify the DA model with a contrast sensitive term as follows:
% $$ r(t) = \frac{\alpha K_{r}\otimes s(t)}{1+\beta (K_{\mu}\otimes{s(t)} + \gamma R(K_{\sigma}\otimes{s(t)} ))} $$
% where $K_{r}$ is a low-pass response filter, $K_{\mu}$ is a low pass gain filter that is sensitive to the mean and is responsible for gain scaling with mean, and $K_{sigma}$ is a high-pass filter that is responsible for gain scaling with contrast. $R$ is the ramp function. 

clear p
p.    A = 157.8750;
p.    B = 4.4062;
p.    C = 2;
p.tau_r = 27.4688;
p.tau_m = 200;
p.tau_s = 20;
p.   s0 = 0;
characteriseModel(@DAModel_contrast,p,data);

prettyFig('fs=14;')


if being_published
	snapnow
	delete(gcf)
end

%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(dataHash(strcat(mfilename,'.m'),Opt))

%%
% This file should be in this commit:
[status,m]=unix('git rev-parse HEAD');
if ~status
	disp(m)
end

t = toc;

%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))

% tag the file as being published 

if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end
