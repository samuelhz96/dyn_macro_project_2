% Project 2, Dynamic Macroeconomics with Numerics
% Hashem Zehi, Samuel (120112285)
% Kotiers, Róza (11945569)
% Polzin, Julian (11948952)
% 18/06/2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  General Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specification of the function which is used to maximize the
% discounted, two-period utility of the agent

% Load parameters:
parameters;
disp(P);



% Create anonymus function:
fun = @(x) P.beta.*(P.alpha.*(P.alpha.*(P.delta.*P.phi).^(-1).*x.^(P.alpha-1)).^(P.alpha./(P.phi-P.alpha)).*x.^(P.alpha-1)+1-P.delta.*(P.alpha.*(P.delta.*P.phi).^(-1).*x.^(P.alpha-1)).^(P.phi./(P.phi-P.alpha)))-1;


% In order to find a suitable starting value for the
% root finding procedure we plot the function defined above
% for a wide range of inputs

% Define vector of consumption values for the plot
c = 1:0.1:100;

% Plot the function defined above:
%plot(c,fun(c));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Root Finding Procedure %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set starting value for iterative method
x0 = 50;


% Call root-finder and print three results
    % 1. x value at which the root is
    % 2. value of the derivative at the root
    % 3. Check for reason of convergence (ideal: 1)
[x,fval,exitflag] = fzero(fun,x0);
x


% Calculate the steady state utilization
U = (P.alpha.*(P.delta.*P.phi).^(-1).*x.^(P.alpha-1)).^(1./(P.phi-P.alpha));
U

% Find the steady-state depreciation
deltaBar = P.delta.*U.^(P.phi);
deltaBar

%Find the steady state output
ybar = (x.*U)^(P.alpha);
ybar
% y in steady state: 3.7471

% Find the steady state consumption
ibar = deltaBar.*x;
ibar
% i in steady state is 0.8993

% Find the steady state investment
cbar = ybar - ibar;
cbar;
% c in steady state is 2.8478

% Verify that the output equals consumption plus investment:
ybar == cbar + ibar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Project 2: Part a: State Space Approach Policy Function %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save the steady state value for capital from the root finding method
kss = x;

% Shortcut for often used term
adpphi = (P.alpha./(P.delta.*P.phi)).^(P.phi./(P.phi-P.alpha));
adpalpha = (P.alpha./(P.delta.*P.phi)).^(P.alpha./(P.phi-P.alpha));

% Steady state consumption as function of steady state capital and parameters
cSS = @(kss) kss.^(P.alpha+(((P.alpha-1).*P.alpha)/(P.phi-P.alpha))).*adpalpha-P.delta.*adpphi.*kss.^(1+(((P.alpha-1).*P.phi)/(P.phi-P.alpha)));

% Steady state derivative of output w.r.t. capital as function of steady state capital and parameters
dySS = @(kss) P.alpha.*adpalpha.*kss.^(P.alpha-1+(((P.alpha-1).*P.alpha)./(P.phi-P.alpha)));

% Steady state depreciation as function of steady state capital and parameters
deltaSS = @(kss) P.delta.*adpphi.*kss.^(((P.alpha-1).*P.phi)/(P.phi-P.alpha));


% Step 1: 
% Partial derivatives of consumption in period t at steady state
dctdkttt = @(kss) 0.*kss;
dctdktt = @(kss) (-kss)./kss;
dctdkt = @(kss) (P.alpha+(((P.alpha-1).*P.alpha)./(P.phi-P.alpha))).*kss.^(P.alpha-1+(((P.alpha-1).*P.alpha)/(P.phi-P.alpha))).*adpalpha+1-P.delta.*adpphi.*(1+(((P.alpha-1).*P.phi)/(P.phi-P.alpha))).*kss.^(((P.alpha-1).*P.phi)/(P.phi-P.alpha));
dctdztt = @(kss) 0.*kss;
dctdzt = @(kss) (1-P.alpha+(((1-P.alpha).*P.alpha)./(P.phi-P.alpha))).*adpalpha.*kss.^(P.alpha+(((P.alpha-1).*P.alpha)/(P.phi-P.alpha)))-P.delta.*adpphi.*(((1-P.alpha).*P.phi)/(P.phi-P.alpha)).*kss.^(1+(((P.alpha-1).*P.phi)/(P.phi-P.alpha)));

% Step 2: 
% Partial derivatives of period t+1 consumption at steady state
dcttdkttt = @(kss) (-kss)./kss;
dcttdktt = @(kss) (P.alpha+(((P.alpha-1).*P.alpha)./(P.phi-P.alpha))).*kss.^(P.alpha-1+(((P.alpha-1).*P.alpha)./(P.phi-P.alpha))).*adpalpha+1-P.delta.*(1+(((P.alpha-1).*P.phi)/(P.phi-P.alpha))).*adpphi.*kss.^(((P.alpha-1).*P.phi)/(P.phi-P.alpha));
dcttdkt = @(kss) 0.*kss;
dcttdztt = @(kss) (1-P.alpha+(((1-P.alpha).*P.alpha)/(P.phi-P.alpha))).*adpalpha.*kss.^(P.alpha+(((P.alpha-1).*P.alpha)/(P.phi-P.alpha)))-P.delta.*adpphi.*(((1-P.alpha).*P.phi)/(P.phi-P.alpha)).*kss.^(1+(((P.alpha-1).*P.phi)/(P.phi-P.alpha)));
dcttdzt = @(kss) P.rho.*(1-P.alpha+(((1-P.alpha).*P.alpha)/(P.phi-P.alpha))).*adpalpha.*kss.^(P.alpha+(((P.alpha-1).*P.alpha)/(P.phi-P.alpha)))-P.rho.*P.delta.*adpphi.*(((1-P.alpha).*P.phi)/(P.phi-P.alpha)).*kss.^(1+(((P.alpha-1).*P.phi)/(P.phi-P.alpha)));

% Step 3:
% Partual derivatives of the f_{k_{t+1}}-delta_{t+1} at steady state
dfdkttt = @(kss) 0.*kss;
dfdktt = @(kss) P.alpha.*(P.alpha-1+(((P.alpha-1).*P.alpha)/(P.phi-P.alpha))).*adpalpha.*kss.^(P.alpha-2+(((P.alpha-1).*P.alpha)/(P.phi-P.alpha))) ;
dfdkt = @(kss) 0.*kss;
dfdztt = @(kss) P.alpha.*(1-P.alpha+(((1-P.alpha).*P.alpha)/(P.phi-P.alpha))).*adpalpha.*kss.^(P.alpha-1+(((P.alpha-1).*P.alpha)/(P.phi-P.alpha)));
dfdzt = @(kss) P.alpha.*P.rho.*(1-P.alpha+(((1-P.alpha).*P.alpha)/(P.phi-P.alpha))).*adpalpha.*kss.^(P.alpha-1+(((P.alpha-1).*P.alpha)/(P.phi-P.alpha)));

% Partial derivative of delta_{t+1} w.r.t. k_t in steady state
ddtdkttt = @(kss) 0.*kss;
ddtdktt = @(kss) (((P.alpha-1).*P.delta.*P.phi)/(P.phi-P.alpha)).*adpphi.*kss.^((((P.alpha-1).*P.phi)/(P.phi-P.alpha))-1);
ddtdkt = @(kss) 0.*kss;
ddtdztt = @(kss) (((1-P.alpha).*P.delta.*P.phi)/(P.phi-P.alpha)).*adpphi.*kss.^((((P.alpha-1).*P.phi)/(P.phi-P.alpha)));
ddtdzt = @(kss) (((1-P.alpha).*P.rho.*P.delta.*P.phi)/(P.phi-P.alpha)).*adpphi.*kss.^(((P.alpha-1).*P.phi)/(P.phi-P.alpha));



% Construct the elements for the matrices of the State Space Approach:
V1 = @(kss) dcttdkt(kss) - P.beta.*(cSS(kss).*(dfdkt(kss)-ddtdkt(kss)) + (dySS(kss)+1-deltaSS(kss)).*dctdkt(kss));
V2 = @(kss) dcttdktt(kss) - P.beta.*(cSS(kss).*(dfdktt(kss)-ddtdktt(kss)) + (dySS(kss)+1-deltaSS(kss)).*dctdktt(kss));
V3 = @(kss) dcttdkttt(kss) - P.beta.*(cSS(kss).*(dfdktt(kss)-ddtdktt(kss)) + (dySS(kss)+1-deltaSS(kss)).*dctdkttt(kss));
V4 = @(kss) dcttdzt(kss) - P.beta.*(cSS(kss).*(dfdzt(kss)-ddtdzt(kss)) + (dySS(kss)+1-deltaSS(kss)).*dctdzt(kss));
V5 = @(kss) dcttdztt(kss) - P.beta.*(cSS(kss).*(dfdztt(kss)-ddtdztt(kss)) + (dySS(kss)+1-deltaSS(kss)).*dctdztt(kss));


A = [V3(kss) V2(kss) V5(kss); 0 1 0; 0 0 1];
B = [0 V1(kss) V4(kss); -1 0 0; 0 0 -P.rho];

PI = -inv(A)*B;
PIabsEigenvalues = abs(eig(PI))
PIeigenvalues = eig(PI,'matrix')

