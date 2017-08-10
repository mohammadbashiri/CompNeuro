clear all; %clear all variables
close all; %close all open figures

t_total = 500; %ms
tstep = 1; %ms
nNeuron = 360;

%variable initialization
Vi = zeros(nNeuron,t_total); %net activation
dSdt = zeros(nNeuron,1);
%Si = zeros(nNeuron,t_total); %Synaptic drive
Fi = zeros(nNeuron,t_total); % Firing rate
extInput = zeros(nNeuron,t_total);
synInput = zeros(nNeuron,t_total);
synInput_V = zeros(nNeuron,t_total);
tonicInhibition = -1.5;
tau = 10;
g = 0.5;

%Wieght matrix generation - gaussian
weightMatGaussian = zeros(nNeuron,nNeuron);
weightMatMexican = zeros(nNeuron,nNeuron);
for i = 1:nNeuron
    for j = 1:nNeuron
        weightMatGaussian(i,j)= weightGenGaussian((i/nNeuron)*2*pi,(j/nNeuron)*2*pi);
        weightMatMexican(i,j)= weightGenMexican((i/nNeuron)*2*pi,(j/nNeuron)*2*pi);
    end
end

weightMatMexican = weightMatGaussian + weightMatMexican;

% weightMat = weightMatMexican;
% weightMat = weightMat - 0.6;

weightMat = weightMatGaussian;

figure;

%plot weight matrix-gaussian
for i=1:10:nNeuron
    hold on;
    plot(weightMat(i,:));
    xlim([0 nNeuron]);
end

%extInput;
inputN = 1;
%extInput(1:50,1) = 11*gausswin(50);
extInput(180:190,1:20) = 10;
extInput(200:210,100:150) = 18;
extInput(100:120,100:150) = 20;
synInput(180,1:10) = 5;

figure;
for t=1:t_total
    
    %computing synaptic input
    synInput_V = synInput;
    for i=1:nNeuron
        Vi(i,t) = tonicInhibition + extInput(i,t) + weightMat(i,:)*synInput_V(:,t);
    end
    
    %computing net activation membrane potential
    %Vi(:,t) = tonicInhibition + extInput(:,t) + synInput(:,t);
    
    %firing rate
    Fi(:,t) = 1./(1+exp(-g*Vi(:,t)));
    %Fi(:,t) = (1+tanh(Vi(:,t)))/2; 
    
    
    dSdt(:,1) = (Fi(:,t)-synInput(:,t))/tau;
    synInput(:,t+1) = synInput(:,t) + dSdt(:,1);
%     
    
%     extInput(inputN+1:inputN+50+1,t+1) = extInput(inputN:inputN+50,t);
%     inputN = inputN + 1;
    
    stem(Fi(1:nNeuron,t));
    xlim([0 nNeuron])
    ylim([0 1.2])
    pause(0.000001);
    t
end
    
%Single Weight Generation function - Gaussian
function weight = weightGenGaussian(theta1,theta2)
    sigma = pi/18;
    if((theta1-theta2)<-pi)
        weight = exp((-(2*pi+theta1-theta2)^2)/(4*sigma^2));
    elseif((theta1-theta2)>pi)
        weight = exp((-(theta1-theta2-2*pi)^2)/(4*sigma^2));
    else
        weight = exp((-(theta1-theta2)^2)/(4*sigma^2));
    end
    weight = (1-weight)*-0.2+weight; %scaling the weight matrix between 1 and - 0.2 - this is without mexican
    weight = weight - 0.5;
    
    %weight = 0.9*weight - 0.5;
    
    %weight = weight*1.5;
end

%Single Weight Generation function - Mexican hat
function weight = weightGenMexican(theta1,theta2)
    sigma = 2*pi/18;
    if((theta1-theta2)<-pi)
        weight = exp((-(2*pi+theta1-theta2)^2)/(4*sigma^2));
    elseif((theta1-theta2)>pi)
        weight = exp((-(theta1-theta2-2*pi)^2)/(4*sigma^2));
    else
        weight = exp((-(theta1-theta2)^2)/(4*sigma^2));
    end
    weight = -weight/2;
end
