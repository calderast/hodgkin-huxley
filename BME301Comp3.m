% Stephanie Crater
% BME 301 Computation 3
% March 22, 2020

clear;

%% Set up variables to investigate runtime vs accuracy
x = linspace(-4,-2,10);
dts = 10.^x; %ms

runs_per_dt = 10; % average these to reduce variability in runtime
runtimes = zeros(1,runs_per_dt); % runtimes to average for each dt 

avg_runtimes = zeros(1,length(dts)); % store avg runtime for each dt
spiketimes = zeros(1,length(dts)); % store spike time for each dt
errors = zeros(1,length(dts)); % store avg error for each dt

for tstep = 1:length(dts)
    for subrun = 1:runs_per_dt % average multiple runs for each dt
    tic;
    %% Define parameters
    Cm = 1.0; %uF/cm^2
    Vrest = -60.0; %mV
    E_Na = 55.17; %mV
    E_K = -72.14; %mV
    E_l = -49.24; %mV
    g_Na = 120.0; %mS/cm^2
    g_K = 36.0; %mS/cm^2
    g_l = 0.3; %mS/cm^2
    
    %% Set options for the simulation
    Istim = 30; %uA/cm^2
    dur = 0.5; %ms
    tmax = 10; %ms
    dt = dts(tstep); %ms
    
    %% Set up time vector and preallocate variables
    t = 0:dt:tmax; %ms
    stim = zeros(1,length(t)); %uA
    Vm = zeros(1,length(t)); %mV
    m = zeros(1,length(t));
    h = zeros(1,length(t));
    n = zeros(1,length(t));
    I_Na = zeros(1,length(t)); %uA/cm^2
    I_K = zeros(1,length(t)); %uA/cm^2
    I_l = zeros(1,length(t)); %uA/cm^2
    
    %% Define functions for gating variables
    am = @(v) 0.1*(-35-v)/(exp((-35-v)/10)-1);
    Bm = @(v) 4*exp((-60-v)/20);
    ah = @(v) 0.07*exp((-60-v)/20);
    Bh = @(v) 1/(exp((-30-v)/10)+1);
    an = @(v) 0.01*(-50-v)/(exp((-50-v)/10)-1);
    Bn = @(v) 0.125*exp((-60-v)/80);
    
    %% Initialize variables
    Vm(1) = Vrest; %mV
    m(1) = am(Vrest)/(am(Vrest)+Bm(Vrest));
    h(1) = ah(Vrest)/(ah(Vrest)+Bh(Vrest));
    n(1) = an(Vrest)/(an(Vrest)+Bn(Vrest));
    
    %% Time loop
    for i = 1:(length(t)-1)
        % Apply stimulus
        if t(i) <= dur
            stim(i) = Istim;
        else
            stim(i) = 0;
        end
        
        % Compute ionic currents
        I_Na(i) = (g_Na*m(i)^3)*h(i)*(Vm(i)-E_Na);
        I_K(i) = (g_K*n(i)^4)*(Vm(i)-E_K);
        I_l(i) = g_l*(Vm(i)-E_l);
        
        % Update state variables
        Vm(i+1) = Vm(i) - dt/Cm*(I_Na(i) + I_K(i) + I_l(i) - stim(i));
        m(i+1) = m(i) + dt*( am(Vm(i))*(1-m(i)) - Bm(Vm(i))*m(i) );
        h(i+1) = h(i) + dt*( ah(Vm(i))*(1-h(i)) - Bh(Vm(i))*h(i) );
        n(i+1) = n(i) + dt*( an(Vm(i))*(1-n(i)) - Bn(Vm(i))*n(i) );
    end
    
    spiketimes(tstep) = t(Vm == max(Vm));
    errors(tstep) = abs(spiketimes(1)-spiketimes(tstep))/spiketimes(1)*100;
    
    runtimes(subrun) = toc;
    end
    % average run time of all trials at a specific dt
    avg_runtimes(tstep) = mean(runtimes);
end

% Plot Vm and Istim
figure(1); clf;
yyaxis left;
plot(t, Vm);
xlabel('Time (ms)');
ylabel('Membrane voltage (mV)');
yyaxis right;
plot(t, stim);
xlabel('Time (ms)');
ylabel('Stimulus Current (uA/cm^2)');
ylim([-10 40]);

% Plot Ik and INa
figure(2); clf;
plot(t, I_K, t, I_Na);
xlabel('Time (ms)');
ylabel('Ion Currents (uA/cm^2)');
legend('I_K', 'I_{Na}');

% PLot m, h and n
figure(3); clf;
plot(t, m, t, h, t, n);
xlabel('Time (ms)');
ylabel('Prob. Gating Variable');
legend('m', 'h', 'n');

% Plot spike time and error vs dt
figure(4); clf;
subplot(2,1,1);
plot(dts, spiketimes);
ylabel('Spike time (ms)');
xlabel('Size of time step (ms)');
hold on;
subplot(2,1,2);
plot(dts, errors);
hold on;
ylabel('Error in spike time (%)');
xlabel('Size of time step (ms)');

% Plot runtime vs dt
figure(5); clf;
plot(dts, avg_runtimes);
ylabel('Program runtime (s)');
xlabel('Size of time step (ms)');

%{
%% Verification exercises
Vm_test = -80.00001:1:60.001;

for V = 1:numel(Vm_test) 
alpham(V) = am(Vm_test(V));
betam(V) = Bm(Vm_test(V));
alphah(V) = ah(Vm_test(V));
betah(V) = Bh(Vm_test(V));
alphan(V) = an(Vm_test(V));
betan(V) = Bn(Vm_test(V));
end

figure(3); clf
subplot(3,1,1);
plot(Vm_test, alpham, Vm_test, betam);
ylabel('kinetics (s^{-1})');
xlabel('Membrane voltage (mV)');
legend('a_m', 'B_m');
xlim([-80 60]);
hold on;
subplot(3,1,2);
plot(Vm_test, alphah, Vm_test, betah);
ylabel('kinetics (s^{-1})');
xlabel('Membrane voltage (mV)');
legend('a_h', 'B_h');
xlim([-80 60]);
subplot(3,1,3);
plot(Vm_test, alphan, Vm_test, betan);
ylabel('kinetics (s^{-1})');
xlabel('Membrane voltage (mV)');
legend('a_n', 'B_n');
xlim([-80 60]);

% Plot Strength vs. duration curve
durations = [0.25 0.375 0.5 0.75 1.0 1.5 2.0 4.0 6.0 8.0 10.0];
strength = [28.14 19.38 14.57 9.99 7.70 5.44 4.35 2.90 2.65 2.62 2.62];
figure(4); clf;
plot(durations, strength);
xlabel('Stimulus duration (ms)');
ylabel('Stimulus strength (uA/cm^2)');
%}
