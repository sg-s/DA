%% data2Filter
% data2Filter is a generalised version of fitFilter2Data, use as:
% K = data2Filter(S,R)
% where S can be a matrix 
% if S is a vector, data2Filter should behave just like fitFilter2Data


function [varargout] = data2Filter(S,R,varargin)



% options and defaults
options.filter_length = 1e3;
options.reg = 0;
options.min_reg = 1e-6;
options.max_reg = 2;
options.reg_max_iter = 10;
options.offset = 200;
options.filter_low_pass = 0;
options.left_trim = 100;
options.right_trim = 100;
options.cross_validate = true;
options.use_cache = true;

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

assert(isvector(R),'2nd argument should be a vector')
assert(length(S) == length(R),'first two arguments should be of equal length')
assert(length(S)>options.filter_length*2,'Data too short')

% orient matrices correctly
if size(S,2)>size(S,1)
	S = S';
end

R = R(:);

% normalise everything
R = R - nanmean(R);
R = R/nanstd(R);

for i = 1:size(S,2)
	S(:,i) = S(:,i) - nanmean(S(:,i));
	S(:,i) = S(:,i)/nanstd(S(:,i));
end

% make the filter struct that we will ultimately return
filterset = struct;
filterset.time = vectorise((1:options.filter_length) - options.offset);
R = circshift(R,options.offset);


% determine verbosity
d = dbstack;
verbosity = 0;
if isempty(find(strcmp('publish',{d.name})))
	verbosity = 1;
end

if ischar(options.reg)
	error('not coded yet. Specify a regularisation factor')
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




only_these_points = options.filter_length+1:length(R);

R = R(only_these_points);

% assemble the matrix equation

% chop up the stimulus into blocks  
reshaped_S = zeros(length(only_these_points),options.filter_length,size(S,2));

for j = 1:size(S,2)
	for i = 1:length(only_these_points)
		reshaped_S(i,:,j) = S(only_these_points(i):-1:only_these_points(i)-options.filter_length+1,j);
	end
end

% combine into a long matrix 
reshaped_S = reshape(reshaped_S,size(reshaped_S,1),size(reshaped_S,3)*size(reshaped_S,2))';

% compute the covariance matrix
C = reshaped_S*reshaped_S';

% regularise over each block
T = sum(reshape(diag(C),length(C)/size(S,2),size(S,2)));

mean_eigen_values = T/options.filter_length;
r = options.reg*mean_eigen_values;


for i = 1:size(S,2)
	a = (i-1)*options.filter_length+1;
	CC = C(a:a+options.filter_length-1,a:a+options.filter_length-1);
	CC = (CC + r(i)*eye(length(CC)))*trace(CC)/(trace(CC) + options.reg*length(CC));
	C(a:a+options.filter_length-1,a:a+options.filter_length-1) = CC;
end

% extract filters by left division 
K = C\(reshaped_S*R);

% split into as many filters as needed
filterset.K = reshape(K,length(K)/size(S,2),size(S,2));


% filter the filters to remove high frequency components, if needed
if options.filter_low_pass > 0
	error('not coded')
	for i = 1:size(S,2)
		filterset.(fn{i}) = filtfilt(ones(options.filter_low_pass,1),options.filter_low_pass,filterset.(fn{i}));
	end
end

% trim filters if need be
filterset.time = filterset.time(options.left_trim+1:end-options.right_trim);
filterset.K = filterset.K(options.left_trim+1:end-options.right_trim,:);

varargout{1} = filterset;




