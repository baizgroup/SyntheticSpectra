function [interpSpec]=dynamicModel(x,aux,t2_fs,interpAx,OutputNoiseShots,SNR)
%won't run until fully updated
t2=t2_fs/1000; %convert from fs to ps
inhomogeneousLW = x(18);
dephasing = x(17); % in ps
LifetimeA = x(15); % in ps.
LifetimeB = x(16); % in ps.
%unpack and convert x variables from experimental values to forms needed
%for model input (cm-1 -> rad/ps, fs -> ps)
%% All peaks
CF=(2*pi)/33.356; %converts from cm-1 to rad/ps
wRot_model=aux.wrot; %Center of exptl. time axis (ps)
delta=-x(1)*CF; %Anharmonicity (cm-1 -> rad/ps)

%% Peak A
tau_a = x(2); %correlation time (ps)
C_a = x(4); %Concentration (can be any units)
u_a = x(7); %Oscillator strength
w_a = (x(10)-aux.wrot).*CF; %frequency, relative to center of the time axis (cm-1 -> rad/ps)
Delta_omega_a = x(12)*CF; %Linewidth (cm-1 -> rad/ps)
%% Peak B
tau_b = x(3); %correlation time (ps)
C_b = x(5); %Concentration (unitless).
u_b = x(8); %Oscillator strength
w_b = (x(11)-aux.wrot).*CF; %frequency, relative to center of the time axis (cm-1 -> rad/ps)
Delta_omega_b = x(13)*CF; %Linewidth (cm-1 -> rad/ps)
%% Crosspeak
kAB=x(6)/(1 + C_a/C_b); %Cross peak rate constant. Comes from kBA/kAB = [A]/[B] and substituting into x(10) = kAB + kBA
kBA=kAB*(C_a/C_b);

%kBA= x(11); %x(11);
Delta_omega_ab=x(14)*CF; %Linewideth (cm-1 -> rad/ps)
%% read dependent variables
%m=In.m; %which loop iteration we're on
t=aux.taxis;

%This function fits a 2D spectrum with 2 peaks, adapted from Hamm and
%Zanni ch. 8

Lambda1 = 1/tau_a;
Lambda2 = 1/tau_b;

%% Define the Kubo Lineshape function
% assume there are some dynamics on the fast modulation limit. Assuming 1ps
% to start
g = @(t) Delta_omega_a.^2/Lambda1.^2.*(exp(-Lambda1.*t)-1+Lambda1.*t) + t/dephasing;
h = @(t) Delta_omega_b.^2/Lambda2.^2.*(exp(-Lambda2.*t)-1+Lambda2.*t) + t/dephasing;
q = @(t) Delta_omega_ab.^2/Lambda1.^2.*(exp(-Lambda1.*t)-1+Lambda1.*t) + t/dephasing;

[T1,T3] = meshgrid(t,t);
%% Define Gamma values,
% these are terms that show up in the differential equations used to model
% the kinetics of vibrational exchange and chemical exchange.

if C_b ~= 0 && C_a ~= 0

    gammaAA = ((kBA*exp(-t2/LifetimeA) + kAB*exp(-(kAB + kBA+1/LifetimeA)*t2))/(kAB + kBA));
    gammaBB = ((kAB*exp(-t2/LifetimeB) + kBA*exp(-(kAB + kBA+1/LifetimeB)*t2))/(kAB + kBA));
    gammaBA = ((kBA*(exp(-t2/LifetimeB) - exp(-(kAB + kBA + 1/LifetimeB)*t2)))/(kAB + kBA));
    gammaAB = ((kAB*(exp(-t2/LifetimeA) - exp(-(kAB + kBA + 1/LifetimeA)*t2)))/(kAB + kBA));

    %% Response functions
    %Peak A
    RA_r = (u_a.^4).*exp(-1i.*w_a.*(-T1+T3)).*gammaAA.*exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*delta.*T3));
    RA_nr = (u_a.^4).*exp(-1i.*w_a.*(T1+T3)).*gammaAA.*exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*delta.*T3));
    %Peak B
    RB_r = (u_b.^4).*exp(-1i*w_b.*(-T1+T3)).*gammaBB.*exp(-h(T1)+h(t2)-h(T3)-h(T1+t2)-h(t2+T3)+h(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*delta.*T3));
    RB_nr = (u_b.^4).*exp(-1i*w_b.*(T1+T3)).*gammaBB.*exp(-h(T1)-h(t2)-h(T3)+h(T1+t2)+h(t2+T3)-h(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*delta.*T3));

    %Cross peak
    RAB_r    = u_a.^2.*u_b.^2.*(exp(-1i*((w_b*T3)-(w_a*T1))).*gammaAB.*exp(-q(T1)).*exp(-q(T3)).*(2-2.*exp(-sqrt(-1)*delta.*T3)));
    RAB_nr   = u_a.^2.*u_b.^2.*(exp(-1i*((w_b*T3)+(w_a*T1))).*gammaAB.*exp(-q(T1)).*exp(-q(T3)).*(2-2.*exp(-sqrt(-1)*delta.*T3)));

    RBA_r    = u_a.^2.*u_b.^2.*(exp(-1i*((w_a*T3)-(w_b*T1))).*gammaBA.*exp(-q(T1)).*exp(-q(T3)).*(2-2.*exp(-sqrt(-1)*delta.*T3)));
    RBA_nr   = u_a.^2.*u_b.^2.*(exp(-1i*((w_a*T3)+(w_b*T1))).*gammaBA.*exp(-q(T1)).*exp(-q(T3)).*(2-2.*exp(-sqrt(-1)*delta.*T3)));

    %the first time points need to be divided by 2 (Sect 9.5.3 in Hamm & Zanni)
    %the reason for this has to do with trapezoidal integration
    RA_r(:,1) = RA_r(:,1)/2;
    RA_r(1,:) = RA_r(1,:)/2;
    RA_nr(:,1) = RA_nr(:,1)/2;
    RA_nr(1,:) = RA_nr(1,:)/2;

    RB_r(:,1) = RB_r(:,1)/2;
    RB_r(1,:) = RB_r(1,:)/2;
    RB_nr(:,1) = RB_nr(:,1)/2;
    RB_nr(1,:) = RB_nr(1,:)/2;

    RAB_r(:,1) = RAB_r(:,1)/2;
    RAB_r(1,:) = RAB_r(1,:)/2;
    RAB_nr(:,1) = RAB_nr(:,1)/2;
    RAB_nr(1,:) = RAB_nr(1,:)/2;

    RBA_r(:,1) = RBA_r(:,1)/2;
    RBA_r(1,:) = RBA_r(1,:)/2;
    RBA_nr(:,1) = RBA_nr(:,1)/2;
    RBA_nr(1,:) = RBA_nr(1,:)/2;

    % add response functions scaled by concentration
    r  = C_a.*RA_r+C_b.*RB_r+C_a.*RAB_r+C_b.*RBA_r; %generate full rephasing spectrum
    nr = C_a.*RA_nr+C_b.*RB_nr+C_a.*RAB_nr+C_b.*RBA_nr; %generate full nonrephasing spectrum

else

    %% Diagonal only
    %Peak A
    RA_r = (u_a.^4).*exp(-1i.*w_a.*(-T1+T3)).*exp(-t2/LifetimeA).*exp(-g(T1)+g(t2)...
        -g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*delta.*T3));
    RA_nr = (u_a.^4).*exp(-1i.*w_a.*(T1+T3)).*exp(-t2/LifetimeA).*exp(-g(T1)-g(t2)...
        -g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*delta.*T3));

    %Peak B
    RB_r = (u_b.^4).*exp(-1i*w_b.*(-T1+T3)).*exp(-t2/LifetimeB).*exp(-h(T1)+h(t2)-...
        h(T3)-h(T1+t2)-h(t2+T3)+h(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*delta.*T3));
    RB_nr = (u_b.^4).*exp(-1i*w_b.*(T1+T3)).*exp(-t2/LifetimeB).*exp(-h(T1)-h(t2)-...
        h(T3)+h(T1+t2)+h(t2+T3)-h(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*delta.*T3));

    %the first time points need to be divided by 2 (Sect 9.5.3 in Hamm & Zanni)
    %the reason for this has to do with trapezoidal integration
    RA_r(:,1) = RA_r(:,1)/2;
    RA_r(1,:) = RA_r(1,:)/2;
    RA_nr(:,1) = RA_nr(:,1)/2;
    RA_nr(1,:) = RA_nr(1,:)/2;

    RB_r(:,1) = RB_r(:,1)/2;
    RB_r(1,:) = RB_r(1,:)/2;
    RB_nr(:,1) = RB_nr(:,1)/2;
    RB_nr(1,:) = RB_nr(1,:)/2;

    r  = C_a.*RA_r+C_b.*RB_r; %generate full rephasing spectrum
    nr = C_a.*RA_nr+C_b.*RB_nr; %generate full nonrephasing spectrum

end

%zero pad the response function
dt = t(2)-t(1);
zp = t(end)+dt:dt:aux.zplength;
zp_t = [t zp]; %zero-padded time axis;

zpmat_r = zeros(length(zp_t));
zpmat_nr = zeros(length(zp_t));

%add the response function to the zero-padded matrix
zpmat_r(1:size(r,1),1:size(r,2)) = r;
zpmat_nr(1:size(nr,1),1:size(nr,2)) = nr;

% do the fft
rephase=flipud(real(fftshift(fft2(zpmat_r))));
nonrephase = rot90(real(fftshift(fft2(zpmat_nr))),2);

modelSpec=rephase+nonrephase;

%normFactor=1/(max(max(modelSpec)));
%Rfinal= modelSpec*normFactor;
% DONT NORMALIZE SPCECTRA
Rfinal= modelSpec;
ax=genFreq(zp_t);
axout = ax + wRot_model;

%zero-pad the noise trajectory along t1
zpNoise = zeros([size(OutputNoiseShots,1) length(zp_t)]);
zpNoise(:,1:size(OutputNoiseShots,2)) = OutputNoiseShots;

NoiseSpec = real(fft(zpNoise,[],2));



% %output the upper quadrant only (faster interpolation);
% Rfinalout= Rfinal(fix(end/2):end,fix(end/2):end);
% axout=ax(fix(end/2):end)+wRot_model;

% output full spectrum (slower interpolation. Only for use with spectra
% covering huge frequency ranges)
interpSpace = (128/(size(interpAx,2)+2));
noiseSpecInterp = interp2(1:128,axout',NoiseSpec', 1:interpSpace:128, interpAx')';

interpSpec_NoNoise = interp2(axout,axout',Rfinal,interpAx, interpAx');
interpSpec_NoNoise = interpSpec_NoNoise./max(max(interpSpec_NoNoise));

%normalize to get appropiate signal to noise
NoiseSpec_norm = (noiseSpecInterp)./std(noiseSpecInterp(:));
NoiseSpec_SNR = NoiseSpec_norm./(SNR.*max(max(abs(interpSpec_NoNoise))));


interpSpec = interpSpec_NoNoise + NoiseSpec_SNR;