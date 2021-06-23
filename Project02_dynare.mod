% 1. Declarations
%--------------------------------------------------------------------------

var k,c,x,U;
varexo e;
parameters beta,alpha,delta,phi,rho;

%--------------------------------------------------------------------------
% 2. Parameter values
%--------------------------------------------------------------------------

beta = 1.03^(-0.25);
alpha = 0.36;
delta = 0.0285;
phi = 1.5;
rho = 0.9;

%--------------------------------------------------------------------------
% 3. Model equations
%--------------------------------------------------------------------------

model;
c(+1) = beta*c*(alpha*exp((1-alpha)*x(+1))*((alpha/(delta*phi))^(1/(phi-alpha))*exp(((1-alpha)/(phi-alpha))*x(+1))*k^((alpha-1)/(phi-alpha)))^(alpha)*k^(alpha-1)+1-delta*((alpha/(delta*phi))^(1/(phi-alpha))*exp(((1-alpha)/(phi-alpha))*x(+1))*k^((alpha-1)/(phi-alpha)))^phi);
k = exp((1-alpha)*x)*(k(-1)*((alpha/(delta*phi))^(1/(phi-alpha))*exp(((1-alpha)/(phi-alpha))*x)*k(-1)^((alpha-1)/(phi-alpha))))^(alpha)+(1-delta*((alpha/(delta*phi))^(1/(phi-alpha))*exp(((1-alpha)/(phi-alpha))*x)*k(-1)^((alpha-1)/(phi-alpha)))^phi)*k(-1)-c;
x = rho*x(-1)+e;
U = (alpha/(phi*delta))^(1/(phi-alpha))*(((1-alpha)*x)/(phi-alpha))*k^((alpha-1)/(phi-alpha));
end;

%--------------------------------------------------------------------------
% 4. Steady states
%--------------------------------------------------------------------------

initval;
c = 2.8478;
k= 60.6232;
x = 0;
U = 0.6;
end;

steady(solve_algo=1);
check;

%--------------------------------------------------------------------------
% 5. Shocks and solution
%--------------------------------------------------------------------------

shocks;
var e;
stderr (0.0005)^(0.5);
end;

stoch_simul(order=1, periods = 1000, irf = 100);
rplot c k U;