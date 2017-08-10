clear all; %clear all variables
close all; %close all open figures

t_total = 10000; %ms
tstep = 1; %ms
nNeuron = 3;
current = [0;2;0]; %zeros(length(inputD),1);
didt = zeros(nNeuron,1);

%weight matrix
v = 0.8;
w = 0;
weightMat(1:nNeuron,1:nNeuron) = v;
for i=1:nNeuron
    weightMat(i,i) = w;
end

weightMat(2,:) = [-0.1 0.2 -0.1];

firing_rate = zeros(nNeuron, t_total);
firing_rate(:,1) = current;%[2 0 1];
tau_s = 100;

weightMat_external = [0.5;0];
external_input = [1 0];

% Introduce the external input
% i think the weights must be the same, however the external input will
% change the dynamics of the system

i1 = animatedline('color', 'r');
i2 = animatedline('color', 'g');
i3 = animatedline('color', 'y');
%i4 = animatedline('color', 'k');
%i5 = animatedline('color', 'm');
axis([1,1000,-4,5]);

for i=1:t_total
    firing_rate(:,i) = current(:,i);
    didt(:,i) = (weightMat * firing_rate(:,i) - current(:,i)) / tau_s;
    %didt(:,i) = (weightMat * firing_rate(:,i) + weightMat_external * external_input(i) - current(:,i)) / tau_s;
    current(:,i+1) = current(:,i) + didt(:,i);
    addpoints(i1,i,current(1,i));
    addpoints(i2,i,current(2,i));
    addpoints(i3,i,current(3,i));
    %addpoints(i4,i,current(4,i));
    %addpoints(i5,i,current(5,i));
    drawnow
    % learn how to plot with handle
end

figure;
hold on;
plot(1:t_total, firing_rate(1,:))
plot(1:t_total, firing_rate(2,:))
plot(1:t_total, firing_rate(3,:))

%activation function