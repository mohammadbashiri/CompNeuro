clear all; %clear all variables
close all; %close all open figures

t_total = 10000; %ms
tstep = 1; %ms
inputD = 2;
outputD = 2;
current = [2;0]; %zeros(length(inputD),1);
didt = zeros(inputD,1);
v = -0.8;
w = 0;
weightMat = [w v;v w];
firing_rate = zeros(inputD, t_total);
%firing_rate(:,1) = [1;1];
tau_s = 100;

firing_rate_extra = zeros(inputD, t_total);

ext_input = 1.0;

for i=1:t_total
    
% this is the Threshold-linear transfer function
%     if current(1,i)>0
%         firing_rate(1,i) = current(1,i);
%     else
%         firing_rate(1,i) = 0;
%     end
%     if current(2,i)>0
%         firing_rate(2,i) = current(2,i);
%     else
%         firing_rate(2,i) = 0;
%     end
    
% this is the quadratic/square root transfer function
    if current(1,i)>1
        firing_rate(1,i) = 2*sqrt(current(1,i)-3/4);
    
    elseif current(1,i)>=0 && current(1,i)<=1
        firing_rate(1,i) = current(1,i)^2;
   
    else
        firing_rate(1,i) = 0;
    end
    
    if current(2,i)>1
        firing_rate(2,i) = 2*sqrt(current(2,i)-3/4);
   
    elseif current(2,i)>=0 && current(2,i)<=1
        firing_rate(2,i) = current(2,i)^2;
    
    else
        firing_rate(2,i) = 0;
    end
    
    firing_rate_extra(2,i) = 0.8*firing_rate(2,i);
    didt(:,i) = (weightMat * firing_rate(:,i) + ext_input - current(:,i)) / tau_s;
    current(:,i+1) = current(:,i) + didt(:,i);
    
%     plot(1:t_total, firing_rate(1,:))
%     plot(1:t_total, firing_rate(2,:))
% 
%     drawnow;
%     pause(1/20);
end
figure;
hold on;
plot(1:t_total, firing_rate(1,:))
plot(1:t_total, firing_rate(2,:))

%activation function

%current = [2;0];
ext_input = zeros(inputD, t_total);
ext_input(1,:) = 0.0;
ext_input(2,:) = 0.0;

for i=1:t_total
    
% this is the Threshold-linear transfer function
%     if current(1,i)>0
%         firing_rate(1,i) = current(1,i);
%     else
%         firing_rate(1,i) = 0;
%     end
%     if current(2,i)>0
%         firing_rate(2,i) = current(2,i);
%     else
%         firing_rate(2,i) = 0;
%     end
    
% this is the quadratic/square root transfer function
    if current(1,i)>1
        firing_rate(1,i) = 2*sqrt(current(1,i)-3/4);
    
    elseif current(1,i)>=0 && current(1,i)<=1
        firing_rate(1,i) = current(1,i)^2;
   
    else
        firing_rate(1,i) = 0;
    end
    
    if current(2,i)>1
        firing_rate(2,i) = 2*sqrt(current(2,i)-3/4);
   
    elseif current(2,i)>=0 && current(2,i)<=1
        firing_rate(2,i) = current(2,i)^2;
    
    else
        firing_rate(2,i) = 0;
    end
    
    firing_rate_extra(2,i) = 0.8*firing_rate(2,i);
    didt(:,i) = (weightMat * firing_rate(:,i) + ext_input(:,i) - current(:,i)) / tau_s;
    current(:,i+1) = current(:,i) + didt(:,i);
    
%     plot(1:t_total, firing_rate(1,:))
%     plot(1:t_total, firing_rate(2,:))
% 
%     drawnow;
%     pause(1/20);
end
figure;
hold on;
plot(1:t_total, firing_rate(1,:))
plot(1:t_total, firing_rate(2,:))