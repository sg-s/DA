% temp file

data.PID = [];
data.voltage = [];


sz = size(PID);
data.PID = zeros(1,sz(3)*sz(1));
for i = 1:sz(2)
	thispid = [squeeze(PID(1,1,:)) ; squeeze(PID(2,1,:)) ; squeeze(PID(3,1,:));  squeeze(PID(4,1,:))];
	thispid = thispid(:);
	data.PID(i,:) = thispid';
	thisorn = [squeeze(ORN(1,1,:)) ; squeeze(ORN(2,1,:)) ; squeeze(ORN(3,1,:));  squeeze(ORN(4,1,:))];
	thisorn = thisorn(:);
	data.voltage(i,:) = thisorn';
end

SamplingRate = 1e4;
ControlParadigm.Outputs = 0*thispid;
ControlParadigm.Name = 'BinaryFlicker';
OutputChannelNames = {'dummy'};

save('ab2A_ethyl_butyrate_fly6_s3.mat','data','ControlParadigm','SamplingRate','OutputChannelNames')