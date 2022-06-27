% This script loads Teliatry data and calculates hemoglobin values and TOI
% values. 

clc
close all
clear

%load NIRS data file
load Hypoxias_MAP_Shout
% load coefficients for hemoglobin and TOI calculation
load e_coef

c_beta = [-1.2132   -0.0728    1.8103    1.1433  -11.5816];

% Calcualte coefficient matrix
wavelength = [730 680 760 850 910]; %Wavelength values provided by Vishnou.
waveidx = round((wavelength - 599.9)/.1);



HBcoef(:,1) = e_coef(waveidx,3);
HBcoef(:,2) = e_coef(waveidx,2);
A = -pinv(HBcoef);

% Calculating hemoglobin concentration values
% W1, ..., W5 are raw wavenlegth values reported by the sensor at each time
% instant. B is the baseline value. If baseline removal is not needed, set
% B=0.
B = 0;
W = nirsData(1:5, :)';
data = (W-B)/16383; % This value (16383) might change. Check with Vishnou.
data = max(data, 0.001)';
HBconc = A*log(data);
O2Hb = HBconc(1,:)';
HHb  = HBconc(2,:)';
tHb  = O2Hb + HHb;

% Calculating TOI. L1, ...,L5 are driving intensities. 
led = double(led);
Amps_coef = [led(1), led(2), led(4), led(5)]/127*16383; % Replace 255*4096 with 127*16383 for Teliatry sensor
L_coefreq = [W(:,1), W(:,2), W(:,4), W(:,5)] ./ W(:,3);

OD_test = -log(L_coefreq./Amps_coef)';
TOI = c_beta*[OD_test; ones(1, size(OD_test,2))];
TOI = abs(TOI')*100;

% Plotting results. 
fs = 20;
time = (1: length(O2Hb))*1/fs;

figure;
plot(time, O2Hb, 'r', time, HHb, 'b', time, tHb, 'g')
legend('O2Hb', 'HHb', 'tHb')
xlabel('Time (s)')

figure;
plot(time, TOI, 'k')
xlabel('Time (s)')
ylabel('TOI')
