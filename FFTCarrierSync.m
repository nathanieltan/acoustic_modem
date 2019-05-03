% Open the file containing the received samples
clear; clf;

Fs = 44100;
Fc = 5000;
period = 100;
sample_offset = period/2;

rx = audioread('rx.wav');
rx = rx(1:length(rx));

rxI = zeros(length(rx),1);
rxQ = zeros(length(rx),1);

% Break received signal into real and imaginary components
for k = 1:length(rx)
    rxI(k) = rx(k) * cos(2*pi*(k/Fs)*Fc);
    rxQ(k) = 1i*rx(k) * -sin(2*pi*(k/Fs)*Fc);
end

% Add real and complex components together
y = rxI + rxQ;
y = lowpass(y, 3000, Fs);

% Trim audio of low-amplitude beginning and end sections
y_start = 1;
y_end = length(y);
keep_looking = true;
abs_mean = abs(mean(real(y)));

for n = 1:length(y)
    if keep_looking && (abs(y(n)) >4.5e2 * abs_mean)
        y_start = n;
        keep_looking = false;
    end
end

keep_looking = true;

for n = length(y):-1:1
    if keep_looking && (abs(y(n)) > 6e2 * abs_mean)
        y_end = n;
        keep_looking = false;
    end
end

% Trimmed signal
y_pkt = y((y_start):(y_end));

mag_h_est = rms(abs(y_pkt));

% Break out phase information using Flat Fading Model
y_est = y_pkt/mag_h_est;
s_time = y_est.^4;
s_freq = fft(s_time);
[spikeVal, fourDelta] = max(s_freq);
delta = -1*(2*pi*(length(s_freq)-1)/length(s_freq))*(fourDelta-1)/((length(s_freq))*4);
theta = -angle(spikeVal)/4;

x_hat = zeros(floor(length(y_est)/period),1);
for k = 1:length(y_est)/period
    x_hat(k) = y_est(k*period - sample_offset)*exp(1i*(delta*((k*period - sample_offset)-1)+theta));
end

% Apply rotation to estimate transmitted values
rotation = exp(1i*pi/4);
x_hat = x_hat.*rotation;

x_adjust = zeros(length(x_hat),1);

% snap estimated values to +- 1 +- j
for k = 1:length(x_hat)
   if real(x_hat(k)) < 0
      x_adjust(k) = x_adjust(k) - 1; 
   elseif real(x_hat(k)) > 0
      x_adjust(k) = x_adjust(k) + 1;
   end
   
   if imag(x_hat(k)) < 0
       x_adjust(k) = x_adjust(k) - 1i;
   elseif imag(x_hat(k)) > 0
       x_adjust(k) = x_adjust(k) + 1i;
   end
       
end

plot(real(x_hat))
title('Real xhat');
figure 
plot(imag(x_hat))
title('Imaginary xhat');
figure
%{
plot(real(y_pkt))
title('Real ypkt');
figure
%}
plot(real(x_hat),imag(x_hat),'.')

title('Constellation');
figure


plot(real(x_adjust));
title('Real Adjust');

audiowrite('output.wav',x_adjust,Fs);
%plot(imag(x_adjust));
%title('Imaginary Adjust');
%figure
