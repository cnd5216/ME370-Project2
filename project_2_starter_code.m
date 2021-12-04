
%% starter code for project 2
clear all; % make sure we don't overwrite other things, so wipe everything clean
close all; % close out prior figures

%% load in data
% this presumes that you are running this m-file in the same folder as the
% file that is loaded
load('mower_deck_base_displacement.mat');
% fs, sampling frequency Hz
% seconds, s, number of seconds of data
% time, time vector in units of s
% base_motion, mm of base displacement at mower deck that holds the mower
% handlebar

%% parameters
wall_thickness=0.0254/16; % m, wall thickness of tubular handlebar (1/16 inch thick)
handlebar_diameter=1*0.0254; % m, outer diameter of tubular handlebar (1 inch diameter)
handlebar_length=40*.0254; % m, length of handlebar (40 inch length)
handlebar_width=18*.0254; % m, width of handlebar at top (18 inch width)
E_handlebar=200e9; % Pa, young's modulus of steel handlebar
rho_handlebar=7200; % kg/m^3, density of steel handlebar
I_handlebar=pi/2*(handlebar_diameter^4/16-(handlebar_diameter-wall_thickness)^4/16); % m^4, moment of inertia of tubular handlebar
A_handlebar=pi/4*(handlebar_diameter^2/4-(handlebar_diameter-wall_thickness)^2/4); % m^2, cross-sectional area of tubular handlebar
mass_handle=handlebar_width*rho_handlebar*A_handlebar; % kg, lumped mass at end of handlebar from the handle system
meq=0.77+mass_handle+2*0.224*rho_handlebar*A_handlebar*handlebar_length; % kg, total dynamic mass of handlebar, 1.7 lb worth of plastic grips and levers, plus metal handlebar pieces, plus the handlebar cantilever dynamic mass component for both sides of the cantilevered handlebar
keq=3*2*E_handlebar*I_handlebar/handlebar_length^3; % N/m, equivalent stiffness of handlebar

%% without absorbers, the handlebar vibrates like a base-excited cantilever
% meq*\ddot{x} + keq*(x - y) = 0
% where y is the base displacement provided in the starter .mat file.
% where x is the displacement of the handlebar mass
% system starts at rest, so no initial conditions
% states are z(1) = x(t) handlebar displacement
% z(2) = \dot{x} handlebar velocity
% state equations are z1dot = z(2)
% z2dot = -keq/meq*(z(1)-base_motion)
% note that base motion is given in mm, so it needs to be converted to m units.
% also because base motion is a predefined vector the
% interp1(time,base_motion,t) code below allows us to use this vector in
% the context of an ode45 simulation that is indexed in an actual time t

[t,y]=ode45(@(t,z)[z(2);-keq/meq*(z(1)-1e-3*interp1(time,base_motion,t))],time,[0;0]);

%% plot untreated handlebar vibration as baseline
figure(1);
clf;
plot(t,1e3*y(:,1),'r');
xlabel('time, s');
ylabel('handlebar displacement, mm');
set(gca,'fontsize',18)

%% fft computation parameters
% the y(:,1) from ode45 simulation output MUST be the handlebar
% displacement state y(:,1) = z(1) = x(t) for this code to work.
% the example code above uses this notation. please copy this 
time_sampled_per_fft=1; % [s] seconds of data acquisition over which fft is evaluated
spacing_cts=round(fs*time_sampled_per_fft); % number of samples to use in FFT to obtain freq_spacing
fft_numbers=floor(fs*seconds/spacing_cts); % number of ffts to compute/loop through
nft=2^nextpow2(spacing_cts); % number of samples next to power of 2 for spacing_cts

%% evaluate frequency response
for ooo=1:2*fft_numbers-1 %  1:2*d.fft_numbers-1 with half-overlap as defined below in trunc
    trunc=(ooo-1)*spacing_cts/2+1:ooo*spacing_cts/2+spacing_cts/2; % define truncation in time : gives half-overlap of windowed averages is the ooo=1:2*d.fft_numbers-1
    yfft=fft(y(trunc,1).*window(@hann,length(trunc)),nft)/(spacing_cts*mean(window(@hann,length(trunc))));
    ch_ft(:,ooo)=2*abs(yfft(1:nft/2+1)); % magnitude of single-sided fourier transform of mower handlebar displacement
end

ch_ft_a=mean(squeeze(ch_ft(1:nft/2+1,:)),2); % average fft of signal
f_ft=fs/2*linspace(0,1,size(ch_ft_a,1))';

%% plot frequency responses
% untreated frequency response of handlebar displacement
figure(2);
clf;
plot(f_ft,1e3*ch_ft_a,'r');
xlabel('frequency, Hz');
ylabel('handlebar displacement amplitude, mm');
set(gca,'fontsize',18,'xscale','log','yscale','log');
xlim([1 1000]);

% untreated frequency response of handlebar acceleration
% in the frequency domain, we can multiply each displacement amplitude by
% omega^2 = (2*pi*frequency)^2 to get the acceleration frequency response
figure(3);
clf;
plot(f_ft,(2*pi*f_ft).^2.*ch_ft_a,'r');
xlabel('frequency, Hz');
ylabel('handlebar acceleration amplitude, m/s^2');
set(gca,'fontsize',18,'xscale','log','yscale','log');
xlim([1 1000]);

%% find RMS accelerations of handlebar
RMSaccel=(sum(((2*pi*f_ft).^2.*ch_ft_a).^2/2)).^(1/2); % m/s^2, RMS of the handlebar displacement frequency response
% we want to reduce the RMSaccel by the application of vibration absorber[s]

%% try one representative handlebar vibration absorber

%%CHANGE

% for a annular rubber bearing, the equivalent radial spring constant is
% k_a = G[shear modulus] * A[cross-sectional area] / L[length
absorber_length=1*.0254; % m, length of absorber piece (4 inch)
absorber_spring_annular_thickness=0.025*.0254; % m, thickness of absorber spring (1/4 inch)
absorber_G=400e3; % Pa, shear modulus of absorber spring
% absorber_A is the cross-sectional area of the annulus of the absorber spring
absorber_A=pi*((handlebar_diameter/2-wall_thickness)^2-(handlebar_diameter/2-wall_thickness-absorber_spring_annular_thickness)^2); % m^2, cross-sectional area of absorber spring
k_a=absorber_G*absorber_A/absorber_length; % N/m, equivalent spring constant of absorber spring
zeta_a=0.07; % damping ratio of absorber spring
absorber_mass_rho=7800; % kg/m^3 density of absorber top mass, steel
absorber_mass_diameter=handlebar_diameter-2*wall_thickness-2*absorber_spring_annular_thickness; % kg, thickness of absorber top mass (1/16 inch thick)
m_a=absorber_mass_rho*pi/4*absorber_mass_diameter^2*absorber_length; % kg, total mass of absorber mass
c_a=zeta_a*2*sqrt(k_a*m_a); % N.s/m damping constant of absorber

tuning_freq=sqrt(k_a/m_a)/2/pi; % Hz tuning frequency 
%%disp(['tuning frequency' num2str(tuning_freq,'10.3f') 'Hz']); %% fix parentheses

%% with a vibration absorber

%CHANGE 

% meq*\ddot{x} + keq*(x - y) + c_a*(\dot{x} - \dot{z}) + k_a*(x - z) = 0
% m_a*\ddot{z} - c_a*(\dot{x} - \dot{z}) - k_a*(x - z) = 0
% where y is the base displacement provided in the starter .mat file.
% where x is the displacement of the handlebar mass
% where z is the displacement of the absorber mass
% system starts at rest, so no initial conditions

[t,y]=ode45(@(t,z)[z(2);-keq/meq*(z(1)-1e-3*interp1(time,base_motion,t))-c_a/meq*(z(2)-z(4))-k_a/meq*(z(1)-z(3));z(4);c_a/m_a*(z(2)-z(4))+k_a/m_a*(z(1)-z(3))],time,[0;0;0;0]);

%% plot treated handlebar vibration with one vibration absorber
figure(1);
hold on
plot(t,1e3*y(:,1),'b');
xlabel('time, s');
ylabel('handlebar displacement, mm');
set(gca,'fontsize',18)
legend('without absorber: baseline','with 1 absorber');

%% evaluate frequency response
for ooo=1:2*fft_numbers-1 %  1:2*d.fft_numbers-1 with half-overlap as defined below in trunc
    trunc=(ooo-1)*spacing_cts/2+1:ooo*spacing_cts/2+spacing_cts/2; % define truncation in time : gives half-overlap of windowed averages is the ooo=1:2*d.fft_numbers-1
    yfft=fft(y(trunc,1).*window(@hann,length(trunc)),nft)/(spacing_cts*mean(window(@hann,length(trunc))));
    ch_ft(:,ooo)=2*abs(yfft(1:nft/2+1)); % magnitude of single-sided fourier transform of mower handlebar displacement
end

ch_ft_a=mean(squeeze(ch_ft(1:nft/2+1,:)),2); % average fft of signal
f_ft=fs/2*linspace(0,1,size(ch_ft_a,1))';

%% plot frequency responses
% untreated frequency response of handlebar displacement
figure(2);
hold on;
plot(f_ft,1e3*ch_ft_a,'b');
xlabel('frequency, Hz');
ylabel('handlebar displacement amplitude, mm');
set(gca,'fontsize',18,'xscale','log','yscale','log');
xlim([1 1000]);

% untreated frequency response of handlebar acceleration
% in the frequency domain, we can multiply each displacement amplitude by
% omega^2 = (2*pi*frequency)^2 to get the acceleration frequency response
figure(3);
hold on
plot(f_ft,(2*pi*f_ft).^2.*ch_ft_a,'b');
xlabel('frequency, Hz');
ylabel('handlebar acceleration amplitude, m/s^2');
set(gca,'fontsize',18,'xscale','log','yscale','log');
xlim([1 1000]);

%% find RMS accelerations of handlebar
RMSaccel_withabsorber=(sum(((2*pi*f_ft).^2.*ch_ft_a).^2/2)).^(1/2); % m/s^2, RMS of the handlebar displacement frequency response
% we want to reduce the RMSaccel by the application of vibration absorber[s]

% display comparison of RMS accelerations
disp(['baseline RMS accel ' num2str(RMSaccel,'%10.3f') ' m/s^2']);
disp(['RMS accel with absorber ' num2str(RMSaccel_withabsorber,'%10.3f') ' m/s^2']);

%% check that absorber displacement does not exceed a physical limit created by the handlebar interior dimensions
max_absorber_disp=max(abs(y(:,3)));
if max_absorber_disp<=handlebar_diameter/2-wall_thickness-absorber_mass_diameter/2
    disp(['No concern: max absorber displacement is less than handlebar constraint']);
else
    disp(['WARNING: max absorber displacement exceeds handlebar constraint']);
end

%%