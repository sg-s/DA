% fitBiophysicalModel
% 
% created by Srinivas Gorur-Shandilya at 4:59 , 18 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


% grab the data
load /code/da/data/MSG_per_neuron
clear d
stim = mean2([MSG_data(1,:).stim]);
resp = mean2([MSG_data(1,:).resp]);
d.stimulus = stim;
d.response = resp;
d.response(1:1e3) = NaN;


p = getModelParameters(@biophysicalModelEuler);
p = fitModel2Data(@biophysicalModelv2,d,'p0',p,'nsteps',200);


