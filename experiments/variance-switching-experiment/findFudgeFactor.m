% findFudgeFactor.m
% finds the best fudge factor so that when we switch variances, the mean remains the same
% 
% created by Srinivas Gorur-Shandilya at 1:25 , 09 December 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

fudge_factor_range = (0.5:0.1:1);

% first, use manual control to make sure you have depleted the headspace
% nicely. 

alldata = struct;

for i = 1:length(fudge_factor_range)
    fudge_factor = fudge_factor_range(i);
    
    % make control paradigm with this
    ControlParadigm = makeVarianceSwitching(fudge_factor);
    
    % now run the experiment
    alldata(i).data = kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',[2 2],'w',10000);
    
end