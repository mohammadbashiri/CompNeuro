function du = intfire_neuronRC(t,u,I_t,param,varargin)
    
    
    
    %if only one I_t is given, assume that I is the same for all neurons
    if length(I_t) == 1
        I_t = ones(size(u)) * I_t;
    end   
    
    R = param(1);
    C = param(2);
    
    du = I_t/C - u/(R*C);

