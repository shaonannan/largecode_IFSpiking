function [LGNt] = temporal_kernel(tau0,tau1,dt,tspan,tdelay)
% [LGNt] = temporal_kernel(0.014,0.056,dt,tfinal);
% [LGNt] = temporal_kernel(0.014,0.056,dt,tfinal,0.01);
if nargin < 2
    error('more parameters needed');
elseif nargin < 3
    dt    = 1e-4;
    tspan = 0.20;   % default spanning time
    tdelay = 0;
elseif nargin < 4
    tspan = 0.20;
    tdelay = 0;
end

tt = dt:dt:tspan;
if tdelay == 0
    LGNt(1,:) = tt.*exp(-tt/tau0)/tau0/tau0 - tt.*exp(-tt/tau1)/tau1/tau1;
    LGNt(2,:) = tt.*exp(-tt/tau0)/tau0/tau0 - tt.*exp(-tt/tau1)/tau1/tau1;
else 
    LGNt(1,:) = tt.*exp(-tt/tau0)/tau0/tau0 - tt.*exp(-tt/tau1)/tau1/tau1;
    tt_d      = (tt-tdelay).*heaviside(tt-tdelay);
    LGNt(2,:) = tt_d.*exp(-tt_d/tau0)/tau0/tau0 - tt_d.*exp(-tt_d/tau1)/tau1/tau1;
end