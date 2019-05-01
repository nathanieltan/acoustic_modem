Fc = 5000;

load handel;
Fs = 8192;
m = sign(y);

audiowrite('test0.wav',m,Fs);

for k = 1:length(x)
    x(k) = m(k) * cos(2*pi*(k/8192)*Fc);
end

audiowrite('test1.wav',x,8192);

y = x;
for k = 1:length(y)
    y(k) = y(k) * cos(2*pi*(k/8192)*Fc);
end

x_est = lowpass(y,3000, 8192);
x_est = sign(x_est);

audiowrite('test2.wav',x_est,8192);