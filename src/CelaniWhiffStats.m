% functions used by Celani et al. to model the distribution of whiff intensities, whiff durations and blank durations in Celani, A., Villermaux, E. & Vergassola, M. Odor Landscapes in Turbulent Environments. Phys. Rev. X 4, 041015–17 (2014).

classdef CelaniWhiffStats < model


	properties
		% define parameters and bounds
		parameter_names = {'A1','A2',  ' B',  'C'};
		lb = 			  [1e-6  1e-6   1e-3  1e-6];
		ub = 			  [1e3   1e3    1e3   1e3 ];

		default_values = [ .001     10    1     1];

		variable_names = {'pS','pD','pB'}; 

	end % end properties 

	methods

		function [m] = evaluate(m)

			x = linspace(min(m.stimulus),max(m.stimulus),1e6);
			m.prediction.pS = (m.parameters.A1./x).*exp(-x./m.parameters.A2);
			m.prediction.pS = m.prediction.pS/sum(m.prediction.pS);
			m.prediction.pS = m.prediction.pS/mean(diff(x));


		end % end evaluate 

		function m = plotDistributions(m,action)
			if ~isfield(m.handles,'plot_fig')
				% this is being called for the first time
				% create a figure
				m.handles.plot_fig = figure('position',[50 250 1500 500],'NumberTitle','off','IntegerHandle','off','Name','Whiff Statistics','CloseRequestFcn',@m.quitManipulateCallback);

				% we need to make only one plot -- 
				
				m.handles.plot_ax = subplot(1,3,1); hold on
				m.handles.plot_ax.XScale = 'log';
				m.handles.plot_ax.YScale = 'log';
				m.handles.plot_ax.XLim = [.01 10];
				xlabel(m.handles.plot_ax,'Whiff intensity')
				ylabel(m.handles.plot_ax,'p(Whiff Intensity)')
				hold(m.handles.plot_ax,'on')
				m.handles.plot_data.handles(1) = plot(m.handles.plot_ax,NaN,NaN,'Color','r');
				prettyFig();
				
			end
				
			if nargin == 2
				if strcmp(action,'update')

					x = linspace(min(m.stimulus),max(m.stimulus),1e6);
					m.evaluate;

					m.handles.plot_data(1).handles(1).XData = x;
					m.handles.plot_data(1).handles(1).YData = m.prediction.pS;
					
					try
						m.handles.plot_ax.YLim = [min(m.prediction.pS) 1];
					catch
					end
				end
			end
		end



 	end % end methods



end % end classdef 