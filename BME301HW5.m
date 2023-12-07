% Stephanie Crater
% BME 301 Homework 5
% February 19, 2020

clear;
%% Set up variables
Vm = linspace(-50, 200, 150);
an = zeros(1,length(Vm));
am = zeros(1,length(Vm));
ah = zeros(1,length(Vm));
Bn = zeros(1,length(Vm));
Bm = zeros(1,length(Vm));
Bh = zeros(1,length(Vm));
n = zeros(1,length(Vm));
m = zeros(1,length(Vm));
h = zeros(1,length(Vm));
dt = .003;
n(1) = .35;
m(1) = .075;
h(1) = .6;

%% Problem 1a: Steady-state values
for i = 1:length(Vm)
    an(i) = 0.01*(10-Vm(i))/(exp((10-Vm(i))/10)-1);
    am(i) = 0.1*(25-Vm(i))/(exp((25-Vm(i))/10)-1);
    ah(i) = 0.07*exp(-Vm(i)/20);
    
    Bn(i) = .125*exp(-Vm(i)/80);
    Bm(i) = 4*exp(-Vm(i)/18);
    Bh(i) = 1/(exp((30-Vm(i))/10)+1);
    
    n(i) = an(i)/(an(i)+Bn(i));
    m(i) = am(i)/(am(i)+Bm(i));
    h(i) = ah(i)/(ah(i)+Bh(i));
end

% Make plots of a and B values
figure(1); clf;
subplot(3,1,1);
plot(Vm, an);
hold on;
plot(Vm, Bn);
legend('a_n', 'B_n');
ylabel('Kinetics (s^{-1})');
title('n');
subplot(3,1,2);
plot(Vm, am);
hold on;
plot(Vm, Bm);
legend('a_m', 'B_m');
ylabel('Kinetics (s^{-1})');
title('m');
subplot(3,1,3);
plot(Vm, ah);
hold on;
plot(Vm, Bh);
xlabel('V (mV)');
legend('a_h', 'B_h');
ylabel('Kinetics (s^{-1})');
title('h');

% Make plots of m,n,h
figure(2); clf;
plot(Vm, n);
hold on;
plot(Vm, m);
plot(Vm, h);
title('Simulation of HH Gating Parameters');
legend('n', 'm', 'h');
xlabel('Vm (mV)');
ylabel('Prob. Gating Variable');

%% Problem 1b: Polyfit
Vrange = linspace(-10,40,50); % shortened voltage range

% get coefficients for approximation
coefs_an = polyfit(Vrange,log(an(40:89)),1);
coefs_am = polyfit(Vrange,log(am(40:89)),1);
coefs_ah = polyfit(Vrange,log(ah(40:89)),1);
coefs_Bn = polyfit(Vrange,log(Bn(40:89)),1);
coefs_Bm = polyfit(Vrange,log(Bm(40:89)),1);
coefs_Bh = polyfit(Vrange,log(Bh(40:89)),1);

% get approximated values
approx_an = exp(coefs_an(2)).*exp(coefs_an(1).*Vrange);
approx_am = exp(coefs_am(2)).*exp(coefs_am(1).*Vrange);
approx_ah = exp(coefs_ah(2)).*exp(coefs_ah(1).*Vrange);
approx_Bn = exp(coefs_Bn(2)).*exp(coefs_Bn(1).*Vrange);
approx_Bm = exp(coefs_Bm(2)).*exp(coefs_Bm(1).*Vrange);
approx_Bh = exp(coefs_Bh(2)).*exp(coefs_Bh(1).*Vrange);

% make plots
figure(3); clf;
subplot(3,2,1);
plot(Vrange, approx_an); hold on;
plot(Vrange, an(40:89), 'k--'); title('a_n');
legend('Boltzman Approx.', 'HH Model');

subplot(3,2,3);
plot(Vrange, approx_am); hold on;
plot(Vrange, am(40:89), 'k--'); title('a_m');
legend('Boltzman Approx.', 'HH Model');

subplot(3,2,5);
plot(Vrange, approx_ah); hold on;
plot(Vrange, ah(40:89), 'k--'); title('a_h');
legend('Boltzman Approx.', 'HH Model');

subplot(3,2,2);
plot(Vrange, approx_Bn); hold on;
plot(Vrange, Bn(40:89), 'k--'); title('B_n');
legend('Boltzman Approx.', 'HH Model');

subplot(3,2,4);
plot(Vrange, approx_Bm); hold on;
plot(Vrange, Bm(40:89), 'k--'); title('B_m');
legend('Boltzman Approx.', 'HH Model');

subplot(3,2,6);
plot(Vrange, approx_Bh); hold on;
plot(Vrange, Bh(40:89), 'k--'); title('B_h');
legend('Boltzman Approx.', 'HH Model');

%% Problem 1c: steady-state m,n,h with Boltzman
n_approx = zeros(1,length(Vrange));
m_approx = zeros(1,length(Vrange));
h_approx = zeros(1,length(Vrange));

for i = 1:50
n_approx(i) = approx_an(i)/(approx_an(i) + approx_Bn(i));
m_approx(i) = approx_am(i)/(approx_am(i) + approx_Bm(i));
h_approx(i) = approx_ah(i)/(approx_ah(i) + approx_Bh(i));
end

figure(4); clf;
plot(Vrange, n(40:89), 'r--');
hold on;
plot(Vrange, m(40:89), 'r:');
plot(Vrange, h(40:89), 'r-');
plot(Vrange, n_approx, 'k--');
plot(Vrange, m_approx, 'k:');
plot(Vrange, h_approx, 'k-');
title('HH and Estimated Gating Parameters');
legend('n (HH)', 'm (HH)', 'h (HH)', 'n (Boltzman)', 'm (Boltzman)', 'h (Boltzman)');
xlabel('Vm (mV)')
ylabel('Prob. Gating Variable');
