% Project 1, Dynamic Macroeconomics with Numerics
% Exercise 2(b)
% Hashem Zehi, Samuel (120112285)
% Kotiers, Róza (11945569)
% Polzin, Julian (11948952)
% 28.05.2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  General Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specification of the function which is used to maximize the
% discounted, two-period utility of the agent

% Load parameters:
parameters;
disp(P)



% Create anonymus function:
fun = @(x) P.beta.*(P.alpha.*(P.alpha.*(P.delta.*P.phi).^(-1).*x.^(P.alpha-1)).^(P.alpha./(P.phi-P.alpha)).*x.^(P.alpha-1)+1-P.delta.*(P.alpha.*(P.delta.*P.phi).^(-1).*x.^(P.alpha-1)).^(P.phi./(P.phi-P.alpha)))-1


% In order to find a suitable starting value for the
% root finding procedure we plot the function defined above
% for a wide range of inputs

% Define vector of consumption values for the plot
c = 1:0.1:100;

% Plot the function defined above:
plot(c,fun(c));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Root Finding Procedure %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set starting value for iterative method
x0 = 50;


% Call root-finder and print three results
    % 1. x value at which the root is
    % 2. value of the derivative at the root
    % 3. Check for reason of convergence (ideal: 1)
[x,fval,exitflag] = fzero(fun,x0)
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
adp = P.alpha./(P.delta.*P.phi);

% Steady state consumption as function of steady state capital and parameters
cSS = @(kss) kss.^(P.alpha+(((P.alpha-1).*P.alpha)/(P.phi-P.alpha))).*(adp).^(P.alpha./(P.phi-P.alpha))-P.delta.*(adp).^(P.phi./(P.phi-P.alpha)).*kss.^(1+(((P.alpha-1).*P.phi)/(P.phi-P.alpha)));

% Steady state output as function of steady state capital and parameters
dySS = @(kss) P.alpha.*adp.^(P.alpha./(P.phi-P.alpha)).*kss.^(P.alpha-1+(((P.alpha-1).*P.alpha)./(P.phi-P.alpha)));

% Steady state depreciation as function of steady state capital and parameters
deltaSS = @(kss) P.delta.*adp.^(P.phi./(P.phi-P.alpha)).*kss.^(((P.alpha-1).*P.phi)/(P.phi-P.alpha));

% Step 1: 
% Partial derivatives of consumption in period t+1 at steady state
dcttdkt = @(kss) 0.*kss         % not actually needed but for completeness
dcttdktt = @(kss) (P.alpha+(((P.alpha-1).*P.alpha)/(P.phi-P.alpha))).*adp.^(P.alpha./(P.phi-P.alpha)).*kss.^(P.alpha-1+(((P.alpha-1).*P.alpha)./(P.phi-P.alpha)))-P.delta.*adp.^(P.phi./(P.phi-P.alpha)).*(1+(((P.alpha-1).*P.phi)/(P.phi-P.alpha))).*kss.^(((P.alpha-1).*P.phi)/(P.phi-P.alpha))
dcttdkttt = @(kss) -kss./kss;
dcttdzt = @(kss) 0.*kss;
dcttdztt = @(kss) (1-P.alpha+(((1-P.alpha).*P.alpha)/(P.phi-P.alpha))).*adp.^(P.alpha./(P.phi-P.alpha)).*kss.^(P.alpha+(((P.alpha-1).*P.alpha)/(P.phi-P.alpha)))-P.delta.*adp.^(P.phi/(P.phi-P.alpha)).*(((1-P.alpha).*P.phi)/(P.phi-P.alpha)).*kss.^(1+(((P.alpha-1).*P.phi)/(P.phi-P.alpha)));

% Step 2: 
% Partial derivatives of consumption in period t at steady state
dctdkt = @(kss) (P.alpha+(((1-P.alpha).*P.alpha)/(P.phi-P.alpha))).*adp.^(P.alpha./(P.phi-P.alpha)).*kss.^(P.alpha-1+(((P.alpha-1).*P.alpha)/(P.phi-P.alpha)))+1-delta.*adp.^(P.phi./(P.phi-P.alpha)).*((((P.alpha-1).*P.phi)/(P.phi-P.alpha)));
dctdktt = @(kss) -kss./kss;
dctdkttt = @(kss) 0.*kss;
%dctdzt = @(kss) 
dctdztt = @(kss) 0.*kss;


% TODO: everything from here on plus one from above!









