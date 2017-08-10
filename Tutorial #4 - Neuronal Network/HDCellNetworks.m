clear all; %clear all variables
close all; %close all open figures

t_total = 1000; %ms
tstep = 1; %ms
nNeuron = 360;
nNetwork = 2;

%variable initialization
Vi = zeros(nNeuron,t_total,nNetwork); %net activation
dSdt = zeros(nNeuron,1,nNetwork);
%Si = zeros(nNeuron,t_total); %Synaptic drive
Fi = zeros(nNeuron,t_total,nNetwork); % Firing rate
extInput = zeros(nNeuron,t_total,nNetwork);
synInput = zeros(nNeuron,t_total,nNetwork);
synInput_V = zeros(nNeuron,t_total,nNetwork);
tonicInhibition = -1.5;
tau = 5;
g = 0.4;

%Wieght matrix generation - gaussian
weightMatGaussian = zeros(nNeuron,nNeuron,nNetwork);
weightMatMexican = zeros(nNeuron,nNeuron,nNetwork);
weightMatplusTheta = zeros(nNeuron,nNeuron,nNetwork);
weightMatminusTheta = zeros(nNeuron,nNeuron,nNetwork);

for i = 1:nNeuron
    for j = 1:nNeuron
        weightMatGaussian(i,j,:)= weightGenGaussian((i/nNeuron)*2*pi,(j/nNeuron)*2*pi);
        weightMatMexican(i,j,:)= weightGenMexican((i/nNeuron)*2*pi,(j/nNeuron)*2*pi);        
    end
end

weightMatMexican = 1.5*weightMatGaussian + weightMatMexican;

% weightMat = weightMatMexican;
% weightMat = weightMat - 0.6;
weightMat = weightMatGaussian;

%shift weight matrix
theta = 20;
for i = 1:nNeuron
    if((i+theta)>nNeuron)
        weightMatminusTheta(i+theta-nNeuron,:,:) = weightMat(i,:,:);
    else
        weightMatminusTheta(i+theta,:,:) = weightMat(i,:,:);
    end
    
    if((i-theta)<=0)
        weightMatplusTheta(nNeuron-abs(i-theta),:,:) = weightMat(i,:,:);
    else
        weightMatplusTheta(i-theta,:,:) = weightMat(i,:,:);
    end
end

%plot weight matrix-gaussian
figure;
for i=1:30:nNeuron
    hold on;
    plot(weightMat(i,:,1));
    %plot(weightMatplusTheta(i,:));
    %plot(weightMatminusTheta(i,:));
    xlim([0 nNeuron]);
end

extInput(170:190,1:20,1) = 5;
extInput(170:190,1:20,2) = 5;

%angVel = 80;
angVel = -100;

angVelMat = randi([-100 100],1,t_total);
for n=1:5
    for i=1:7-n:t_total-10
        angVelMat(i:i+10) = sum(angVelMat(i:i+10))/10;
    end
end

% figure;
% plot(angVelMat);

%lesion
percentage = 80;
lr = randi([1 nNeuron],1,nNeuron*(percentage/100));
rr = randi([1 nNeuron],1,nNeuron*(percentage/100));

figure;
for t=1:t_total
    
    for i=2:nNeuron-1
        Vi(i,t,1) = tonicInhibition + extInput(i,t,1) + weightMatminusTheta(i,:,1)*synInput(:,t,1)...
        + weightMatplusTheta(i,:,1)*synInput(:,t,2);
    
%         Vi(i,t,1) = tonicInhibition + extInput(i,t,1) + weightMatminusTheta(i,:,1)*synInput(:,t,1); % right movement
%         Vi(i,t,1) = tonicInhibition + extInput(i,t,1) + weightMatplusTheta(i,:,1)*synInput(:,t,1); % left movement
    
        Vi(i,t,2) = tonicInhibition + extInput(i,t,2) + weightMatplusTheta(i,:,2)*synInput(:,t,2)...
             + weightMatminusTheta(i,:,2)*synInput(:,t,1);
    end
    
    %computing net activation membrane potential
    %Vi(:,t) = tonicInhibition + extInput(:,t) + synInput(:,t);
    
    %firing rate
    %angVel = angVel+200/t_total;
%     lVel = 1 + (angVel/100);
%     rVel = 1 - (angVel/100);
    
    lVel = 1 + (angVelMat(t)/100);
    rVel = 1 - (angVelMat(t)/100);
    
    Fi(:,t,:) = 1./(1+exp(-g*Vi(:,t,:)));
    Fi(:,t,1) = Fi(:,t,1)*lVel;
    Fi(:,t,2) = Fi(:,t,2)*rVel;
    %Fi(:,t,:) = (1+tanh(Vi(:,t,:)))/2; 
    
    %lesion
    Fi(rr,t,1) = 0;
    Fi(lr,t,2) = 0;
    
    Fi(1,t,:) = 0;
    Fi(nNeuron,t,:) = 0;
    
    dSdt(:,1,:) = (Fi(:,t,:)-synInput(:,t,:))/tau;
    synInput(:,t+1,:) = synInput(:,t,:) + dSdt(:,1,:);
    
    subplot(2,1,1)
    stem(Fi(1:nNeuron,t,1));
    xlim([0 nNeuron])
    ylim([0 2])
    subplot(2,1,2)
    stem(Fi(1:nNeuron,t,2));
    xlim([0 nNeuron])
    ylim([0 2])
    pause(0.01);
    t
    
    angVel = angVel + 200/t_total; 
end
figure;
subplot(3,1,1);
image(Fi(:,:,1),'CDataMapping','scaled');
subplot(3,1,2);
image(Fi(:,:,2),'CDataMapping','scaled');
subplot(3,1,3);
plot(angVelMat);

%Single Weight Generation function - Gaussian
function weight = weightGenGaussian(theta1,theta2)
    sigma = pi/18;
    if((theta1-theta2)<-pi)
        weight = exp((-(2*pi+theta1-theta2)^2)/(2*sigma^2));
    elseif((theta1-theta2)>pi)
        weight = exp((-(theta1-theta2-2*pi)^2)/(2*sigma^2));
    else
        weight = exp((-(theta1-theta2)^2)/(2*sigma^2));
    end
    
    weight = (1-weight)*-0.2+weight; %scaling the weight matrix between 1 and - 0.2 - this is without mexican
%   weight = weight - 0.4;  
%   %weight = 0.9*weight - 0.5;
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