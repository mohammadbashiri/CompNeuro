%Head direction cell network (e.g. Boucheny et al., 2005), (Song and Wang,
%2005)
%This version treats the population as a whole in the simulation (only one
%Weight matrix)

dt = 0.1;
nTimepoints = 10000;
t = (0:nTimepoints-1)*dt;

v_In = 0.2;
tau = 1;

%synaptic time constant

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
       
E = [ones(numNeurons, 1); -1*ones(numNeurons,1); zeros(numNeurons,1)];
G = [g_Inh*ones(numNeurons, 1); g_Inh*ones(numNeurons,1); g_Ex*ones(numNeurons,1)];



%input to the network (~turning velocity)

v = zeros(1, nTimepoints);
v(1:nTimepoints*0.9) = -1;
v(1:nTimepoints*0.5) = +1;
v(1:nTimepoints*0.1) = 0;
v = v * v_In;


% starting activity
s_DTNl = 0.2*ww;
s_DTNr = s_DTNl;
s_LMN = s_DTNl;
%
f_DTNl = 0.2*ww;
f_DTNr = f_DTNl;
f_LMN = f_DTNl;

s = [s_DTNl; s_DTNr; s_LMN];
f = [f_DTNl; f_DTNr; f_LMN];

figure
hline(1)=line('XData',x,'YData',f_DTNl, 'Color', 'r');
hline(2)=line('XData',x,'YData',f_DTNr, 'Color', 'b');
hline(3)=line('XData',x,'YData',f_LMN, 'Color', 'k');


for n = 1:nTimepoints
    
    vv = G + Weights * s + E*v(n);
    f = max(vv, 0);
    ds = (f - s)/tau;
    
    s = s + ds*dt;
    
    
%     % DTN_l
%     v_DTNl = g_Inh + WE*s_LMN + wi * s_DTNr + v(n);
%     f_DTNl = max(v_DTNl, 0); %Activation function
%     ds_DTNl = (f_DTNl - s_DTNl)/tau;
%     
%     % DTN_r
%     v_DTNr = g_Inh + WE*s_LMN + wi * s_DTNl - v(n);
%     f_DTNr = max(v_DTNr, 0); %Activation function
%     ds_DTNr = (f_DTNr - s_DTNr)/tau;
%      
%     v_LMN = g_Ex + WIL * s_DTNl + WIR * s_DTNr;
%     f_LMN = max(v_LMN, 0); %Activation function
%     ds_LMN = (f_LMN - s_LMN)/tau;
%     
    s_DTNl = s_DTNl + ds_DTNl*dt;
    s_DTNr = s_DTNr + ds_DTNr*dt;
    s_LMN = s_LMN + ds_LMN*dt;
    
    set(hline(1),'YData',f(1:101))
    set(hline(2),'YData',f(102:202))
    set(hline(3),'YData',f(203:303))
    drawnow()
    
    %pause(0.001)
    
    
    %pop(n)=atan2(mean(sin(x/50*pi).*f_LMN'),mean(cos(x/50*pi).*f_LMN'))*50/pi;
end

%figure, plot(pop)