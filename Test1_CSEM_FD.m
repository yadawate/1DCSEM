
%--------------------------------------------------------------------------
clc;
clear all;
close all;
%--------------------------------------------------------------------------
% Adjustable model and data parameters:
%--------------------------------------------------------------------------
global k h srln n_int n_layer omega mu_h eps_h zz offset Rcdpt Srdpt
%       Air    Ocean    Seafloor--->
z = [-1e60 0 300 1300 1400 2400];   % Layer top depths (the first value is not used)
rho = [1e12 0.3125 1 2 1 1];   % Layer resistivities  (ohm-m)
mu0 = pi * 4e-7;
eps0 = 8.8541878176e-12;
Rcdpt = 290;
Srdpt = 270;
f   = 0.25;                       % Frequency (Hz)
zTx = Srdpt;                     % Depth of transmitter (m)
r   = linspace(100,8000,48);  % Ranges to receivers   (m) 
zRx = Rcdpt*ones(size(r));      % Depth of receivers   (m) 
sig = 1./rho;
zz = Rcdpt - Srdpt;
sigma = sig';
mu_h = mu0 * ones(length(sig),1);
eps = eps0 * ones(length(sig),1);
omega = 2 * pi * f;
eps_h = eps + (1j * sigma)/ omega;
h = [4000 300 1000 100 1000 4000];
k = (1+1j)*sqrt(omega*mu0*sigma/2);
%k = omega*sqrt(mu_h.*eps_h);
srln = 2;
n_layer = length(mu_h);
n_int = n_layer - 1;

for rh = 1:length(r)
    offset = r(rh);
    fun1 = @(krho) (krho./kzf(krho,k(srln))).*besselj(0,krho*offset).*(exp(1j*sqrt(k(srln)^2 - krho.^2)*abs(zz))+loseth_te(krho));
    fun2 = @(krho) (1./kzf(krho,k(srln))).*besselj(1,krho*offset).*(exp(1j*sqrt(k(srln)^2 - krho.^2)*abs(zz))+ loseth_te(krho));
    fun3 = @(krho) (krho.*kzf(krho,k(srln))).*besselj(0,krho*offset).*(exp(1j*sqrt(k(srln)^2 - krho.^2)*abs(zz))+loseth_tm(krho));
    fun4 = @(krho) (kzf(krho,k(srln))).*besselj(1,krho*offset).*(exp(1j*sqrt(k(srln)^2 - krho.^2)*abs(zz))+loseth_tm(krho));
    
    IA0TE = omega * mu0 * integral(fun1,0,Inf,'AbsTol',1e-17,'RelTol',1e-17);
    IA1TE = omega * mu0 * integral(fun2,0,Inf,'AbsTol',1e-17,'RelTol',1e-17);
    IA0TM = (omega *mu0/(k(srln).^2)) * integral(fun3,0,Inf,'AbsTol',1e-17,'RelTol',1e-17);
    IA1TM = (omega *mu0/(k(srln).^2)) * integral(fun4,0,Inf,'AbsTol',1e-17,'RelTol',1e-17);

    Erho(rh) = (-1/(4*pi)) * (IA0TM + (IA1TE - IA1TM)/offset);
    Ebet(rh) = (-1 / (4 * pi)) * (-IA0TE + (1 / offset) * (IA1TE - IA1TM));
    
end

figure;
semilogy(r,abs(Erho))
xlabel('Range (m)')
ylabel('Amplitude (V/Am^2)')




