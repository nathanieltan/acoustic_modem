clear;

Fc = 5000;

N = 1000;

% make N random bits of values +- 1
bits = sign(randn(N,1));

Symbol_period = 100;

% create a generic pulse of unit height
% with width equal to symbol period
pulse = ones(Symbol_period, 1);
% Pulse amplitude modulation generation
% Symbol_period is width of pulse and period

% spread out the values in "bits" by Symbol_period
% first create a vector to store the values
x = zeros(Symbol_period*length(bits),1);

% assign every Symbol_period-th sample to equal a value from bits
x(1:Symbol_period:end) = bits;

% now convolve the single generic pulse with the spread-out bits
x_tx = conv(pulse, x);

x_tx = [zeros(1000,1);x_tx;zeros(1000,1)];


Fs = 44100;



audiowrite('test0.wav',bits,Fs);

audiowrite('premod.wav',x_tx,Fs);

for k = 1:length(x_tx)
    x_tx(k) = x_tx(k) * cos(2*pi*(k/Fs)*Fc);
end

audiowrite('test1.wav',x_tx,Fs);

y = x_tx;
for k = 1:length(y)
    y(k) = y(k) * cos(2*pi*(k/Fs)*Fc);
end

x_est = lowpass(y,3000, Fs);
x_est = sign(x_est);

audiowrite('test2.wav',x_est,Fs*Symbol_period);