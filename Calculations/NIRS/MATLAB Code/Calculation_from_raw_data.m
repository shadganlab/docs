load e_coef
load('TOIcoefs2018', 'c_beta');

wavelength = [670 730 810 850 950]; % Replace these values with [730 680 760 850 910] for Teliatry sensor
waveidx = round((wavelength - 599.9)/.1);
HBcoef(:,1) = e_coef(waveidx,3);
HBcoef(:,2) = e_coef(waveidx,2);
A = -pinv(HBcoef);

% Calculating hemoglobin concentration values
% W1, ..., W5 are raw wavenlegth values reported by the sensor at each time
% instant. B is the baseline value. If baseline removal is not needed, set
% B=0.
data = ([W1 W2 W3 W4 W5]-B)/4096; % Replace 4096 with 16383 for Teliatry sensor
data = max(data, 0.001)';
HBconc = A*log(data);
O2Hb = HBconc(1,:)';
HHb  = HBconc(2,:)';
tHb  = O2Hb + HHb;


% Calculating TOI. L1, ...,L5 are driving intensities. 
Amps_coef = [L1, L2, L4, L5]/255*4095; % Replace 255*4096 with 127*16383 for Teliatry sensor
L_coefreq = [W1, W2, W4, W5] ./ W3;

OD_test = -log(L_coefreq./Amps_coef)';
TOI = c_beta*[OD_test; ones(1, size(OD_test,2))];
TOI = abs(TOI')*100;