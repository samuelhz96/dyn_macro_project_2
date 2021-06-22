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
parameters1;
% Load alternative parameters (comment out unused)
%parameters2
%parameters3
%parameters4

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
adp = (P.alpha./(P.delta.*P.phi)).^(1./(P.phi-P.alpha));
adpphi = (P.alpha./(P.delta.*P.phi)).^(P.phi./(P.phi-P.alpha));
adpalpha = (P.alpha./(P.delta.*P.phi)).^(P.alpha./(P.phi-P.alpha));

% Steady state consumption as function of steady state capital and parameters
cSS = @(kss) kss.^(P.alpha+(((P.alpha-1).*P.alpha)/(P.phi-P.alpha))).*adpalpha-P.delta.*adpphi.*kss.^(1+(((P.alpha-1).*P.phi)/(P.phi-P.alpha)));

% Steady state derivative of output w.r.t. capital as function of steady state capital and parameters
dySS = @(kss) P.alpha.*adpalpha.*kss.^(P.alpha-1+(((P.alpha-1).*P.alpha)./(P.phi-P.alpha)));

% Steady state depreciation as function of steady state capital and parameters
deltaSS = @(kss) P.delta.*adpphi.*kss.^(((P.alpha-1).*P.phi)/(P.phi-P.alpha));

% Steady state capital utiliuzation
USS = @(kss) (P.alpha./(P.delta.*P.phi)).^(1./(P.phi-P.alpha)).*kss.^((P.alpha-1)./(P.phi-P.alpha));

% Step 1: DF_1(*)
DF1 = [
    -P.beta.*(P.alpha.*USS(kss).^(P.alpha).*kss.^(P.alpha-1)+1-P.delta.*USS(kss).^P.phi), 0, 0, 0;
    0, 0, 0, 0;
    -1, P.delta.*P.phi.*USS(kss).^(P.phi-1).*kss+P.alpha.*kss.^(P.alpha).*USS(kss).^(P.alpha-1), (1-P.delta.*USS(kss).^P.phi)+P.alpha.*kss.^(P.alpha-1).*USS(kss).^(P.alpha), (1-P.alpha).*(kss.*USS(kss)).^(P.alpha);
    0, 0, 0, -P.rho];

% Step 2: DF_2(*)
DF2 = [1, -P.beta.*cSS(kss).*(P.alpha.^2.*(kss.*USS(kss)).^(P.alpha-1)+1-P.delta.*P.phi.*USS(kss).^(P.phi-1)), -P.beta.*cSS(kss).*P.alpha.*(P.alpha-1).*USS(kss).^(P.alpha).*kss.^(P.alpha-2), -P.beta.*cSS(kss).*P.alpha.*(1-P.alpha).*kss.^(P.alpha-1);
0, 1, -adp.*((P.alpha-1)./(P.phi-P.alpha)).*kss.^(((P.alpha-1)./(P.phi-P.alpha))-1), -adp.*((1-P.alpha)./(P.phi-P.alpha)).*kss.^((P.alpha-1)./(P.phi-P.alpha));
0, 0, -1, 0;
0, 0, 0, 1];

% Step 3: Build Jacobian
Jac = -(inv(DF2))*DF1;

[V,D] = eig(Jac);

eigvals         = diag(D);              % eigenvalues from D as a vector 
idx_stable      = find(abs(eigvals)<1); % Index of stable eigenvalues 
eigvals_stable  = eigvals(idx_stable)   % Retrieves stable eigenvalues fromv vector
eigvecs_stable  = V(:,idx_stable); 
    
ukx_eigvecs = eigvecs_stable(2:4,:);
c_eigvecs = eigvecs_stable(1,:);
   
inv_ukx_eigvecs= inv(ukx_eigvecs);
   
const1= c_eigvecs*inv_ukx_eigvecs(:,1);
const2= c_eigvecs*inv_ukx_eigvecs(:,2);
const3= c_eigvecs*inv_ukx_eigvecs(:,3);
   
G1= Jac(3,1).*const2+Jac(3,3)
G2= Jac(3,1).*const3+Jac(3,4)
G0 = kss

    



