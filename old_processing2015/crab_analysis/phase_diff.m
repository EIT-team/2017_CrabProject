function phase_lag = phase_diff(x,y)
%
% take the FFT
X=fft(x);
Y=fft(y);

% Determine the max value and max point.
% This is where the sinusoidal
% is located. 
[mag_x idx_x] = max(abs(X));
[mag_y idx_y] = max(abs(Y));

% determine the phase difference
% at the maximum point.
px = angle(X(idx_x));
py = angle(Y(idx_y));
phase_lag = py - px;