% meant to be run by GaussianConstructor.m
% 
% created by Srinivas Gorur-Shandilya at 4:18 , 02 March 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [ControlParadigm] = MakeGaussianFlicker(mfc_min,mfc_max,T,dt,tc)
clear ControlParadigm
ControlParadigm.Outputs= zeros(3,T/dt);
ControlParadigm.Name = 'Flicker';
n = floor(T/tc);
stream = RandStream.getGlobalStream;
load rand_state
stream.State = savedState;

nsteps= T/tc;
noise = rand(nsteps,1);
% set to a random number every tc

a = 1;
z = a + tc/dt;
for j = 1:nsteps
    
    ControlParadigm.Outputs(2,a:z) = noise(j)*(mfc_max-mfc_min) + mfc_min;
    

    % increment
    a = z;
    z = min([a + tc/dt T/dt]);
end
    
ControlParadigm.Outputs(2,end) = mean(ControlParadigm.Outputs(2,:)); % don't leave the ORN with some high or low value
ControlParadigm.Outputs(2,ControlParadigm.Outputs(2,:) > 5) = 5; % clip for sanity
ControlParadigm.Outputs(2,ControlParadigm.Outputs(2,:) < 5/200) = 5/200; % clip for sanity

% add the main air
ControlParadigm.Outputs(3,:) = 1;
ControlParadigm.Name = 'all_off';

% make a all zero
ControlParadigm(2).Outputs = zeros(3,1e4);