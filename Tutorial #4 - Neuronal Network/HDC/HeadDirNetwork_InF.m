%Head direction cell network (e.g. Boucheny et al., 2005), (Song and Wang,
%2005)
%This version treats the population as a whole in the simulation (only one
%Weight matrix) and simulates using Integrate and Fire neurons

dt = 0.1;
nTimepoints = 10000;
t = (0:nTimepoints-1)*dt;

v_In = 0.2;
tau = 1;%synaptic time constant

%Number of neurons per population (DTNl, DTNr, LMN)
numNeurons = 101;
shift = round(numNeurons/10);

%create kernel
sd = 0.2;
x = linspace(-1,1,numNeurons);
ww=exp(-x.^2/sd^2)';
ww = circshift(ww, floor(numNeurons/2));

wi = -0.5;

%synaptic connectivity matrix
for i = 1:numNeurons
    W(i,:) = circshift(ww,i);
end

WE=0.2*W;
WI=wi*eye(size(W));%-1.5;
WIL=-5*circshift(W, floor(numNeurons/2)-shift+1);
WIR=-5*circshift(W, floor(numNeurons/2)+shift);

g_Inh = -1.0;
g_Ex = 1.0;


Weights = [zeros(size(W)), WI, WE;...
           WI, zeros(size(W)), WE;...
           WIL, WIR, zeros(size(W))];
       
%input to the network (~turning velocity)

v = zeros(1, nTimepoints);
v(1:nTimepoints*0.9) = -1;
v(1:nTimepoints*0.5) = +1;
v(1:nTimepoints*0.1) = 0;
v = v * v_In;

I = [repmat(v, numNeurons, 1);
    repmat(-v, numNeurons, 1);
    repmat(0*v, numNeurons, 1)];

%I = I + randn(size(I));





u0 = rand(size(Weights,1),1);


[spikes, u] = SolveNeuronODE_EulerRC(t,Weights,u0,I,'firing_thresh',1,'R',Inf,'C',1e-10,'tau_syn',1.05,'syn_prop_delay',2e-3,'psc_shape','exp');