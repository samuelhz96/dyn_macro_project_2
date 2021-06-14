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



(x.*U).^(P.alpha)-((P.beta.*cbar).*(P.alpha.*U.*(x.*U).^(P.alpha-1)+1-P.delta.*U.^(P.phi)))

