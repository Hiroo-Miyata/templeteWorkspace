function [ eout, ec] = ecgRemovalFilter(emg,ecg, emgrest, ecgrest,fs)

if length(emgrest) ~= length(ecgrest)
    error('emgrest and ecgrest data-length must be the same');
end

[ ns, nc ] = size(emgrest);
p = round(0.05*fs); % parameter p = nr. of samples over 0.05 s 

for c = 1:nc % for each emg channel
    U = zeros(ns,p); % matrix of delayed inputs
    U(:,1) = ecgrest;
    for i = 2:p
        U(i:ns,i) = ecgrest(1:ns-i+1);
    end % for i
    h(:,c) = U\emgrest(:,c);  % least square fit of impulse response h
    
end % for c

[ ns, nc ] = size(emg);
ec = zeros(ns,nc);
for c = 1:nc
    ec1 = conv(ecg,h(:,c));% apply convolution
    ec(:,c) = ec1(1:ns);
end
% correction applied and final HP filter
eout = emg-ec;