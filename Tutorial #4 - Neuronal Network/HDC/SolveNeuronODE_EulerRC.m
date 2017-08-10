function [spikes, u] = SolveNeuronODE_EulerRC(time, W_syn,u_0,I,varargin)


%default parameters
S.nrndiffeq = @intfire_neuronRC;
S.R = 1e7; %Ohms
S.C = 1e-8; %Farad
S.firing_thresh = 1; %Arbitrary
S.tau_syn = 0.010; %seconds
S.syn_prop_delay = 0; %seconds
S.psc_shape = 'delta';
S.V_ext = zeros(size(time));
S.minpotential = -Inf;


%custom parameters
n_additional_arguments = numel(varargin);
%change default properties
for n = 1:n_additional_arguments
    if ischar(varargin{n})
        S.(varargin{n}) = varargin{n+1};
        varargin{n+1} = [];
        if ~isfield(S,varargin{n})
            warning(['No default property ' varargin{n}])
        end
    end
end






%simulation parameters

t_stop = max(time);
timestep = mean(diff(time));
if any(abs(diff(time) - timestep) > 10*eps(max(time)))
    error('Time is not sampled with constant timestep!')
end

nrnparams = [S.R, S.C];
tau_m = S.R*S.C;




if any(size(W_syn) ~= length(u_0))
    error('W_syn dimensions and number of starting values do not match!')
end





%initialize u variable
u = zeros(length(u_0),length(time));
u(:,1) = u_0;
%initialize spikes variable
spikes = zeros(size(u,1),length(time));

%postsynaptic current pulse
%(http://icwww.epfl.ch/~gerstner/SPNM/node26.html)

tau_syn = S.tau_syn / timestep; %give tau_syn in timestep units

if strcmp(S.psc_shape,'exp')
    alpha = exp(-(0:timestep:5*S.tau_syn)/S.tau_syn); %alpha has an eponentially decaying shape (postsynaptic current pulse)
elseif strcmp(S.psc_shape,'delta') || strcmp(S.psc_shape,'pulse')
    alpha = 1;
end

alpha = ((alpha / sum(alpha)) * S.C * S.firing_thresh) / timestep; %normalize for "charge conservation"

I = I * S.C;

delay_ts = round(S.syn_prop_delay/timestep);

%still need to fix the units in general: what units are tau and R given in,
%what units are used down in the loop!


pct = 5;
prog = pct/100*length(time);
%%
for t = 1:length(time)-1
    
    %solve membrane potential for the next time step with one euler step
    u(:,t+1) = u(:,t) + timestep * S.nrndiffeq(t,u(:,t),I(:,t),nrnparams,S.V_ext(t));
    
    u(:,t+1) = max(u(:,t+1), S.minpotential);
    
    %find all neurons that fire an action potential
    firing = u(:,t+1)-S.V_ext(t+1) > S.firing_thresh;
    
    %reset potentials and propagate spikes
    if any(firing)
        u(firing,t) = 3*S.firing_thresh; %just for illustration
        %reset neuron membrane potential
        u(firing,t+1) = u(firing,t+1) - S.firing_thresh; %substract threshold instead of setting to 0 so nothing gets lost.
        %what could get lost would depend (for a non-delta-kernel alpha) on
        %the timestep (which we dont want)
        
        spikes(firing,t+1) = 1;    %log spikes
        
        %propagate spikes (with delay) (by modifying the upcoming input
        %currents)
        %I(:,t+1) = I(:,t+1) + W_syn' * firing; %simple Dirac spikes
        
        I_spikes = zeros(size(I));
        synapse_activity = W_syn' * firing;
        synaptic_current = synapse_activity * alpha;
        I_spikes(:,t+1+delay_ts:t+length(alpha)+delay_ts) = synaptic_current;
        I_spikes = I_spikes(1:size(I,1),1:size(I,2)); %truncate everything that is outside the simulation time range
        I = I + I_spikes;
        
        if length(time) > 1e4 && t > prog
            disp([num2str(pct) '% done'])
            pct = pct + 5;
            prog = pct/100*length(time);
        end
        
    end
end
spikes = logical(spikes);

%%

figure('Name','Voltage traces')
subplot(2,1,1)
plot(time,u')
title('Membrane potential')
subplot(2,1,2)
plot(time,I)
title('Membrane current')
%legend('Membrane Potential','Syn Current')
%figure
%plot(I_spikes')
%title('I_{spikes}')

figure('Name','Spikes')
for u = 1:size(spikes,1)
    for s = 1:size(spikes,2)
        if spikes(u,s)
            line([s*timestep s*timestep], [u-1 u],'Color','k')
        end
    end
end
set(gca,'XLim',[min(time) max(time)])
set(gca,'YLim',[0 length(u_0)])


