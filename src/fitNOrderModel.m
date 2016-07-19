%% fitNOrderModel.m
% fits a first order model:
% 0 = K0 x s(t) + (1 + K1 x s(t)) r(t)
% where x is a convolution
% to data specified by s(t) and r(t)
% and K0 and K1 are filters that are fit non-parametrically 


function [varargout] = fitNOrderModel(S,R,varargin)



% options and defaults
options.filter_length = 1e3;
options.reg = 0;
options.min_reg = 1e-6;
options.max_reg = 2;
options.reg_max_iter = 10;
options.offset = 200;
options.filter_low_pass = 10;
options.left_trim = 100;
options.right_trim = 100;
options.cross_validate = true;
options.use_cache = true;
options.model_order = 1;

if nargout && ~nargin 
	varargout{1} = options;
    return
end

% validate and accept options
if iseven(length(varargin))
	for ii = 1:2:length(varargin)-1
	temp = varargin{ii};
    if ischar(temp)
    	if ~any(find(strcmp(temp,fieldnames(options))))
    		disp(['Unknown option: ' temp])
    		disp('The allowed options are:')
    		disp(fieldnames(options))
    		error('UNKNOWN OPTION')
    	else
    		options.(temp) = varargin{ii+1};
    	end
    end
end
elseif isstruct(varargin{1})
	% should be OK...
	options = varargin{1};
else
	error('Inputs need to be name value pairs')
end

assert(isvector(S),'First argument must be a vector')
assert(isvector(R),'2nd argument should be a vector')
assert(length(S) == length(R),'first two arguments should be of equal length')
assert(length(S)>options.filter_length*2,'Data too short')

S = S(:);
% make the filter struct that we will ultimately return
filterset = struct;
filterset.time = vectorise((1:options.filter_length) - options.offset);
S = circshift(S,-options.offset);
R = R(:);

% determine verbosity
d = dbstack;
verbosity = 0;
if isempty(find(strcmp('publish',{d.name})))
	verbosity = 1;
end

if ischar(options.reg)

	if options.cross_validate
		if verbosity
			disp('find best regularization by cross validation')
		end
		z = floor(length(S)/2);
		S_test = S(1:z);
		R_test = R(1:z);
		S = S(z+1:end);
		R = R(z+1:end);
	else
		if verbosity
			disp('find best regularization by self-validation')
		end
		S_test = S;
		R_test = R;
	end


	r2 = NaN(options.reg_max_iter,1);

	% check the cache. 
	temp = options;
	temp.stimulus = S;
	temp.response = R;
	hash = dataHash(temp);
	if isempty(cache(hash))

		% round 1
		if verbosity
			disp('Round 1...')
		end
		reg_vec = logspace(log10(options.min_reg),log10(options.max_reg),options.reg_max_iter);

		for i = 1:options.reg_max_iter
			if verbosity
				textbar(i,options.reg_max_iter)
			end
			options.reg = reg_vec(i);
			temp = fitNOrderModel(S,R,options);
			switch options.model_order
			case 1
				R_pred = firstOrderModel(S_test,temp);
			otherwise
				error('not coded 117')
			end
			
			r2(i) = rsquare(R_pred,R_test);
		end
		[~,idx] = max(r2);
		options.min_reg = reg_vec(idx)/2;
		options.max_reg = reg_vec(idx)*2;

		temp3 = reg_vec(:);
		temp4 = r2(:);

		% round 2
		if verbosity
			disp('Round 2...')
		end
		reg_vec = logspace(log10(options.min_reg),log10(options.max_reg),options.reg_max_iter);

		for i = 1:options.reg_max_iter
			if verbosity
				textbar(i,options.reg_max_iter)
			end
			options.reg = reg_vec(i);
			temp = fitNOrderModel(S,R,options);
			switch options.model_order
			case 1
				R_pred = firstOrderModel(S_test,temp);
			otherwise
				error('not coded 117')
			end
			r2(i) = rsquare(R_pred,R_test);
		end

		diagnostics.reg_vec = [temp3; reg_vec(:)];
		diagnostics.r2 = [temp4; r2(:)];

		% save to cache
		cache(hash,diagnostics)


	else
		if verbosity
			disp('Using cached data...')
		end
		diagnostics = cache(hash);

	end

	options.reg = diagnostics.reg_vec(diagnostics.r2 == max(diagnostics.r2));

	if nargout == 2
		varargout{2} = diagnostics;
	end


end



assert(options.reg >= 0,'regularisation must not be negative')



if options.reg > 0
	S = S + options.reg*randn(length(S),1);
end


only_these_points = options.filter_length+1:length(R);

R = R(only_these_points);

% assemble the matrix equation
switch options.model_order
case 1
	% chop up the stimulus into blocks  
	stim = zeros(options.filter_length,length(only_these_points));

	for i = 1:length(only_these_points)
		stim(:,i) = S(only_these_points(i):-1:only_these_points(i)-options.filter_length+1);
	end

	% also multiply each row by the response and save that too
	stim_r = stim;
	for i = 1:size(stim_r,1)
		stim_r(i,:) = stim_r(i,:).*R';
	end

	% combine these two
	shat = [stim; stim_r];

	C = shat*shat';
	K = C\(shat*R);

	filterset.K0 = K(1:length(K)/2);
	filterset.K1 = -K(length(K)/2+1:end);
end



% filter the filters to remove high frequency components, if needed
fn = fieldnames(filterset);
fn(strcmp(fn,'time')) = [];

if options.filter_low_pass > 0
	for i = 1:length(fn)
		filterset.(fn{i}) = filtfilt(ones(options.filter_low_pass,1),options.filter_low_pass,filterset.(fn{i}));
	end
end

% trim filters if need be
fn = fieldnames(filterset);
for i = 1:length(fn)
	filterset.(fn{i}) = filterset.(fn{i})(options.left_trim+1:end-options.right_trim);
end

varargout{1} = filterset;




