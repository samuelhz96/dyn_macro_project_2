function [dc,dk,dz] = lineartransition(P,dk0,dz0,T)
    % computes the trajectory of (c,k,z) implied by the linearized dynamics
    % starting from (dk0,dz0)=(k0-kss,z0-zss) 
    %
    % inputs:  P   parameters [structure]
    %          dk0 initial StSt-deviation capital stock
    %          dz0 initial StSt-deviation TFP level
    %          T   time horizon
    % outputs: dc  time path of deviations from consumption [T+1 vector]
    %          dk  time path of deviations from capital [T+1 vector]
    %          dz  time path of deviations from TFP [T+1 vector]

    % solution path satisfies dx(t) = a1*e1^t*v1 + a2*e2^t*v2
    %                         dx(0) = (dc0,dk0,dz0)
    
    % Step 1: determine a1 and a2 such that (compare to slide 12)
    %            dk0 = a1*v1(2) + a2*v2(2)
    %            dz0 = a1*v1(3) + a2*v2(3)
    V = [P.v1(2:3) P.v2(2:3)];
    a = V \ [dk0;dz0];               % [dk0;dz0] = a*V
    
    % Step 2: calculate time path (compare to slide 9)
    dx = a(1) * P.v1 * P.e1.^(0:T) + a(2) * P.v2 * P.e2.^(0:T);
    dc = dx(1,:);
    dk = dx(2,:);
    dz = dx(3,:);
    
end
    
    
    