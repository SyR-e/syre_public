% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.


function [f,mag] = FFTAnalysis(y,Fw)
% Performs the FFT analysis

n = length(y);
T = 1/Fw;                               % rotation period
Fs = n/T;                               % sampling frequency
Ts = 1/Fs;                              % sampling period
f = Fs*(0:(n/2))/n;                     % Real frequency
t = (0:n-1)*Ts;                          % Real time vector

yFFT = fft(y);

% Calculates amplitudes
mag = 2*abs(yFFT/n);
mag(1) = mag(1)/2;
mag = mag(1:n/2+1);

% % %  debug below
% n = length(B_Sta_r(:,jj));
% Fs = n*Fw;                               % sampling frequency
% Ts = 1/Fs;                            % sampling period
% f = Fs*(0:n/2)/n;                     % Real frequency
% t = (0:n-1)*Ts*1e3;                   % Time vector [ms]
% 
% figure(1)
% % plot(t,B_Sta_r(:,jj),'r*')
% plot(t,B_Sta_r(:,jj));
% xlabel('One period/ms')
% ylabel('Flux density/T')
% grid on
% % % 
% % 
% yFFT = fft(B_Sta_r(:,jj));
% 
% mag = 2*abs(yFFT/n);
% mag(1) = mag(1)/2;
% mag = mag(1:n/2+1);
% figure(3)
% bar(f,mag);
% xlabel('Frequency/Hz')
% ylabel('Flux density/T')
% xlim([0,2000])
