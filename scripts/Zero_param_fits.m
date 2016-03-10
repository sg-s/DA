% Zero_param_fits.m
% 
% created by Srinivas Gorur-Shandilya at 1:52 , 10 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.



pHeader;

%% Zero Parameter Fits to Data
% In this document we attempt to fit a general DA-like model to data. The model looks like this:
% 
% $$ r(t)=\alpha\left(K_{r}\otimes s(t)\right)\left(\frac{1}{1+\beta_{\mu}K_{\mu}\otimes s(t)}\right)\left(\frac{1}{1+\beta_{\sigma}\left|K_{\sigma}\otimes s(t)\right|}\right) $$
% 

%%
% where the first term controls gain by the mean of the stimulus and the second term controls gain by the standard deviation of the stimulus. Here, we attempt to directly find the values of these parameters from our various experiments, so we can plug them in and see if they can account for responses to naturalistic stimuli. 




prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


