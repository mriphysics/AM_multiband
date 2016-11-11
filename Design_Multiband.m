clear all
close all;
% Code to design either conventional or Amplitude-modulated (AM) only
% Mulitband RF pulses using 
% 1) Phase-optimzation (Wong, ISMRM proceedings 2012)
% 2) Time-shifting (Auerbach et al MRM 2013)
% 3) Root-flipping (Sharma et al MRM 2015)

% Requires CVX (http://cvxr.com/cvx/) and Pauly's RF tools
% (http://rsl.stanford.edu/research/software.html).

gamma_mT = 267522.1; % radians/mT/s
b1max = 20*1e-3; % <--- peak B1 
mb = 6; % <-- Number of slices
tb = 4; % <-- Time-bandwidth product
bs = 10; % <-- Slice separation, in units of slice thicknesses
slthick = 2*1e-3;
AM_only = 0;    %<--- set to 0 for unconstrained design. 1 for AM design.

% Load in single-band refocusing pulses for matched-excitation
load('SB_SLR_cvxdesign_flip180_matchedexcitation.mat');

rfsb = pulse(tb/2).rf;
Nt = 2048;
rfsb = length(rfsb)/Nt*interp1(linspace(0,1,length(rfsb)),rfsb,linspace(0,1,Nt));

FOV = 0.3;
Nz = 4000;
z = linspace(-FOV/2,FOV/2,Nz)';
pos = [0*z 0*z z];
% Design Phase-optimized waveform
rfmb_po = Phaseopt_fn(rfsb,mb,tb,bs,AM_only);

dt_po = max(abs(rfmb_po))/gamma_mT/b1max;

BW = tb/(length(rfmb_po)*dt_po);
Gsel = 2*pi*BW/(gamma_mT*slthick);


figure(1);
Gz = Gsel*ones(length(rfmb_po),1);
G = [0*Gz 0*Gz Gz];

[~,~,~,~,~,b180_po] = blochsim_CK(rfmb_po/gamma_mT/dt_po,G,pos,ones(Nz,1),zeros(Nz,1),'dt',dt_po);
Mxy_po = b180_po(:,end).^2;

% Design Time-shifted waveform

frac_shift = 0.5;
rfmb_ts = Timeshift_fn(rfsb,mb,tb,bs,frac_shift,AM_only);

dt_ts = max(abs(rfmb_ts))/gamma_mT/b1max;

T_ts = (1 + frac_shift)*length(rfmb_ts)*dt_ts;%<-- not sure this is right..

BW = tb/(length(rfsb)*dt_ts);
Gsel = 2*pi*BW/(gamma_mT*slthick);

figure(1);
Gz = Gsel*ones(length(rfmb_ts),1);
G = [0*Gz 0*Gz Gz];

[~,~,~,~,~,b180_ts] = blochsim_CK(rfmb_ts/gamma_mT/dt_ts,G,pos,ones(Nz,1),zeros(Nz,1),'dt',dt_ts);
Mxy_ts = b180_ts(:,end).^2;

% Design Root-flipped waveform
[rfmb_rf90,rfmb_rf180,tbsc]= rootflip_fn(512,mb,tb,bs,AM_only);

dt_rf = max(abs(rfmb_rf180))/(gamma_mT*b1max);
BW = tbsc/(length(rfmb_rf180)*dt_rf);
Gsel = 2*pi*BW/(gamma_mT*slthick);        

Gz = Gsel*ones(length(rfmb_rf180),1);
G = [0*Gz 0*Gz Gz];

[~,~,~,~,~,b180_rf] = blochsim_CK(rfmb_rf180/gamma_mT/dt_rf,G,pos,ones(Nz,1),zeros(Nz,1),'dt',dt_rf);
Mxy_rf = b180_rf(:,end).^2;

% Plot RF waveforms
figure(1);
subplot(311); plot(t_po,rfmb_po_mT );
title('Phase-optimizing');
xlabel('Time [s]');ylabel('B_1 [mT]')
subplot(312); plot(t_ts,rfmb_ts_mT );
title('Time-shifting');
xlabel('Time [s]');ylabel('B_1 [mT]')
subplot(313); plot(t_rf,rfmb_rf_mT );
title('Root-flipping');
xlabel('Time [s]');ylabel('B_1 [mT]')
% plot(abs(rfmb_ts));
% plot(abs(rfmb_rf180));

% Plot refocusing profile |Mxy|
figure(2);
plot(z,abs(Mxy_po));
hold on;grid on;
plot(z,abs(Mxy_ts));
plot(z,abs(Mxy_rf));
legend('Phase-Optimzing','Time-Shifting','Root-Flipping');

t_po = (0:length(rfmb_po)-1)*dt_po;
rfmb_po_mT = abs(rfmb_po)/gamma_mT/dt_po;

t_ts = (0:length(rfmb_ts)-1)*dt_ts;
rfmb_ts_mT = abs(rfmb_ts)/gamma_mT/dt_ts;

t_rf = (0:length(rfmb_rf180)-1)*dt_rf;
rfmb_rf_mT = abs(rfmb_rf180)/gamma_mT/dt_rf;

