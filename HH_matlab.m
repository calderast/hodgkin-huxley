%% Hodgkin-Huxley simulation
%{ 
Some code to simulate an axon's response to an applied stimulus current 
based on the Hodgkin-Huxley equations as described in the 1952 paper 
(cited below). Allows for visualization of how the gating variables 
m, h, and n influence ionic currents and membrane voltage, how the gating 
variables change based on their respective rate constants a and B, and how 
those rate constants are dependent on the membrane voltage.

Provides code to simulate an action potential using the moden convention 
(resting potential at -65mV and positive ions into the cell = voltage gets 
more positive) because that's what we're used to. I replicate the graphs 
from the Hodgkin-Huxley paper using both the modern and the original 
Hodgkin-Huxley convention (resting potential at 0mV and sign of voltage is 
flipped)

In "The Annotated Hodgkin and Huxley" (linked below), the modern versions
of the graphs should match what is simulated here using the modern
convention.

HODGKIN AL, HUXLEY AF. A quantitative description of membrane current and 
its application to conduction and excitation in nerve. 
J Physiol. 1952 Aug;117(4):500-44. doi: 10.1113/jphysiol.1952.sp004764. 
PMID: 12991237; PMCID: PMC1392413.

The Annotated Hodgkin and Huxley: A Reader's Guide
https://press.princeton.edu/books/paperback/9780691220635/the-annotated-hodgkin-and-huxley
%}

%% Define constants
clear;

% constants from Table 3 in the paper
Cm = 1.0; % uF/cm^2
g_Na = 120.0; % mS/cm^2
g_K = 36.0; % mS/cm^2
g_l = 0.3; % mS/cm^2

% set convention = "modern" or "original"/"HH"
convention = "modern";

switch convention

    case "modern"
        %% Define functions for gating variables (modern convention)
        %{
Note that these functions and equilibrium potentials differ from the ones 
in the original Hodgkin-Huxley paper:
The sign of the voltage has been flipped to be consistent with the modern 
convention and the voltage is shifted down by 60mV so the resting potential
of the axon is at -60mV instead of 0mV. I would cite my source for the
Nernst potentials but I forgot where I got them from, sorry
        %}

        Vrest = -65.0; % mV
        E_Na = 55.17; % mV
        E_K = -72.14; % mV
        E_l = -49.24; % mV

        an = @(v) 0.01*(-50-v)/(exp((-50-v)/10)-1);
        Bn = @(v) 0.125*exp((-60-v)/80);
        am = @(v) 0.1*(-35-v)/(exp((-35-v)/10)-1);
        Bm = @(v) 4*exp((-60-v)/18);
        ah = @(v) 0.07*exp((-60-v)/20);
        Bh = @(v) 1/(exp((-30-v)/10)+1);

        Istim = 30; % uA/cm^2

    case {"original", "HH"}
        %% Define functions for gating variables (original HH convention)
        %{
Equilibrium potentials and equations as written in the original
Hodgkin-Huxley paper. The constants come from Table 3.
In "The Annotated Hodgkin and Huxley", equations for the rate constants 
are on page 222 and Table 3 is on page 224. 
        %}

        % Right now the simulation part doesn't work with original HH parameters
        Vrest = 0; % mV
        E_Na = -115; % mV
        E_K = 12; % mV
        E_l = -10.613; % mV

        an = @(v) 0.01*(v+10)/(exp((v+10)/10)-1);
        Bn = @(v) 0.125*exp(v/80);
        am = @(v) 0.1*(v+25)/(exp((v+25)/10)-1);
        Bm = @(v) 4*exp(v/18);
        ah = @(v) 0.07*exp(v/20);
        Bh = @(v) 1/(exp(v/10)+1);

        Istim = -80; % uA/cm^2
end

%% Set options for the simulation
% Stimulus current and duration for when we simulate an action potential

%Istim = -90; % uA/cm^2
dur = 0.5; % ms
tmax = 10; % ms
dt = 10^-3; % ms
t = 0:dt:tmax; % ms

% Preallocate for speed
stim = zeros(1,length(t)); % uA
Vm = zeros(1,length(t)); % mV
m = zeros(1,length(t));
h = zeros(1,length(t));
n = zeros(1,length(t));
I_Na = zeros(1,length(t)); % uA/cm^2
I_K = zeros(1,length(t)); % uA/cm^2
I_l = zeros(1,length(t)); % uA/cm^2

% Initialize variables
Vm(1) = Vrest; % mV
m(1) = am(Vrest)/(am(Vrest)+Bm(Vrest));
h(1) = ah(Vrest)/(ah(Vrest)+Bh(Vrest));
n(1) = an(Vrest)/(an(Vrest)+Bn(Vrest));

%% Simulate the axon's response to a stimulus current

for i = 1:(length(t)-1)
    % Apply stimulus current for specified duration
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
    m(i+1) = m(i) + dt*(am(Vm(i))*(1-m(i)) - Bm(Vm(i))*m(i));
    h(i+1) = h(i) + dt*(ah(Vm(i))*(1-h(i)) - Bh(Vm(i))*h(i));
    n(i+1) = n(i) + dt*(an(Vm(i))*(1-n(i)) - Bn(Vm(i))*n(i));
end

%% Make plots of a simulated action potential

% Plot simulated membrane voltage as a result of the applied stimulus current
figure(1); clf;
yyaxis left;
plot(t, Vm);
xlabel('Time (ms)');
ylabel('Membrane voltage (mV)');
yyaxis right;
plot(t, stim);
xlabel('Time (ms)');
ylabel('Stimulus Current (uA/cm^2)');
%ylim([-10 40]);

% Plot simulated potassium current (I_K) and sodium current (I_Na)
figure(2); clf;
plot(t, I_K, t, I_Na);
xlabel('Time (ms)');
ylabel('Ion Currents (uA/cm^2)');
legend('I_K', 'I_{Na}');

% Plot m, h and n over the course of the action potential
figure(3); clf;
plot(t, m, t, h, t, n);
xlabel('Time (ms)');
ylabel('Prob. Gating Variable');
legend('m', 'h', 'n');

%% Recreate graphs from the Hodgkin-Huxley paper
Vm_test = -80.00001:1:60.001;

alpha_m = zeros(1,length(Vm_test));
beta_m = zeros(1,length(Vm_test));
alpha_h = zeros(1,length(Vm_test));
beta_h = zeros(1,length(Vm_test));
alpha_n = zeros(1,length(Vm_test));
beta_n = zeros(1,length(Vm_test));
m_ = zeros(1,length(Vm_test));
h_ = zeros(1,length(Vm_test));
n_ = zeros(1,length(Vm_test));

for V = 1:numel(Vm_test)
    % Calculate rate constants for each membrane voltage
    alpha_m(V) = am(Vm_test(V));
    beta_m(V) = Bm(Vm_test(V));
    alpha_h(V) = ah(Vm_test(V));
    beta_h(V) = Bh(Vm_test(V));
    alpha_n(V) = an(Vm_test(V));
    beta_n(V) = Bn(Vm_test(V));

    m_(V) = alpha_m(V)/(alpha_m(V)+beta_m(V));
    h_(V) = alpha_h(V)/(alpha_h(V)+beta_h(V));
    n_(V) = alpha_n(V)/(alpha_n(V)+beta_n(V));
end

figure(4); clf
subplot(3,1,1);
plot(Vm_test, alpha_m, Vm_test, beta_m);
ylabel('kinetics (s^{-1})');
xlabel('Membrane voltage (mV)');
legend('a_m', 'B_m');
xlim([-80 60]);
hold on;
subplot(3,1,2);
plot(Vm_test, alpha_h, Vm_test, beta_h);
ylabel('kinetics (s^{-1})');
xlabel('Membrane voltage (mV)');
legend('a_h', 'B_h');
xlim([-80 60]);
subplot(3,1,3);
plot(Vm_test, alpha_n, Vm_test, beta_n);
ylabel('kinetics (s^{-1})');
xlabel('Membrane voltage (mV)');
legend('a_n', 'B_n');
xlim([-80 60]);

figure(5); clf;
subplot(3,1,1);
plot(Vm_test, m_);
ylabel('Prob. Gating Variable');
xlabel('Membrane voltage (mV)');
legend('m');
xlim([-80 60]);
hold on;
subplot(3,1,2);
plot(Vm_test, h_);
ylabel('Prob. Gating Variable');
xlabel('Membrane voltage (mV)');
legend('h');
xlim([-80 60]);
subplot(3,1,3);
plot(Vm_test, n_);
ylabel('Prob. Gating Variable');
xlabel('Membrane voltage (mV)');
legend('n');
xlim([-80 60]);