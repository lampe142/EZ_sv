% Model: RBC model with Epstein Zin Preferences
%   
% Author: Dario Caldara
%
% Purpose: Prepared for the course: "Tools for Nonlinear DSGE models 
% Last Update: 14.06.2010
%
%----------------------------------------------------------------
% 0. Housekeeping
%----------------------------------------------------------------

close all

%----------------------------------------------------------------
% 1. Endogenous variables
%----------------------------------------------------------------
var 

// Allocation variables 
k y c l i

// Utility variables
u v ev 

// Input prices
r w 

rf

// Bond yield variables

// Shocks
z sigma_zt
;
%----------------------------------------------------------------
% 2. Exogenous variables
%----------------------------------------------------------------

varexo e e_sigma;

%----------------------------------------------------------------
% 3. Parameters
%----------------------------------------------------------------

parameters 

// Utility function
nu beta gamma psi theta cte
// Technology 
alpha delta lambda sigma rho eta
;

%----------------------------------------------------------------
% 4. Calibration
%----------------------------------------------------------------

// Utility 
beta = 0.991;
nu =  0.3622;
// gamma = 5; 
gamma = 40; 
psi = 0.5; 
theta = (1-gamma)/(1-1/psi);

// Technology
alpha = 0.3;
delta = 0.0196;
lambda = 0.95;
// sigma = 0.007;
sigma = 0.021;
rho = 0.9;
eta = 0.1;
%----------------------------------------------------------------
% 5. Steady State
%----------------------------------------------------------------

A = ((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1));
B = (nu*(1-alpha)*A^(alpha))/((1-nu)*(A^(alpha-1)-delta));

l_ss = B/(A+B);
k_ss  = B*(1-l_ss);
c_ss = k_ss^alpha*l_ss^(1-alpha)-delta*k_ss;
w_ss = (1-alpha)*k_ss^alpha*l_ss^(-alpha);
r_ss  = 1+alpha*k_ss^(alpha-1)*l_ss^(1-alpha)-delta;
cte = (c_ss^nu*(1-l_ss)^(1-nu))^(-1);

// u_ss = c_ss^nu*(1-l_ss)^(1-nu)
u_ss = cte*c_ss^nu*(1-l_ss)^(1-nu)
 
v_ss = u_ss

ev_ss = v_ss^(1-gamma);

%----------------------------------------------------------------
% 6. Model
%----------------------------------------------------------------

//model(use_dll); 
model;

    // 1. Utitlity
    u = cte*c^nu*(1-l)^(1-nu);

    // 2. Expected Value Function
    ev = v(+1)^(1-gamma);

    // 3. Value Function
    v = ((1-beta)*u^((1-gamma)/theta)+ beta*(ev^(1/theta)))^(theta/(1-gamma));

    // 4. Static Leisure-Consumption
    (1-nu)/nu*c/(1-l)= w;

    // 6. Euler Equation for Capital
    beta*(u(+1)/u)^((1-gamma)/theta)*(c/c(+1))*(v(+1)^(1-gamma)/ev)^(1-1/theta)*r(+1) = 1;
    //beta*(u(+1)/u)^((1-gamma)/theta)*(c/c(+1))*(v(+1)^(1-gamma)/v(+1)^(1-gamma))^(1-1/theta)*r(+1) = 1;
    //beta*(u(+1)/u)^((1-gamma)/theta)*(c/c(+1))*(1)^(1-1/theta)*r(+1) = 1;
    //beta*(u(+1))^((1-gamma)/theta)*(1/c(+1))*(v(+1)^(1-gamma)/ev)^(1-1/theta)*r(+1) = u^((1-gamma)/theta)*(1/c);


    // 8. Production Function
    y = exp(z)*k(-1)^alpha*l^(1-alpha);

    // 9. Return on Risky Asset
    r = 1+alpha*y/k(-1) - delta;

    // 10. Wage
    w = (1-alpha)*y/l;

    // 11. Resource Constraint
    c + i = y;

    // 12. Law of Motion for Capital
    k = (1-delta)*k(-1) + i;

    // 13. Law of Motion for Productivity
    z = lambda*z(-1) + sigma_zt*e;

    // 14. Law of Motion for Volatility
    log(sigma_zt) = (1-rho)*log(sigma) + rho*log(sigma_zt(-1)) + eta*e_sigma;

    // 15. Risk-free Rate
    beta*(u(+1)/u)^((1-gamma)/theta)*(c/c(+1))*(v(+1)^(1-gamma)/ev)^(1-1/theta)*rf = 1;

end;

%----------------------------------------------------------------
% 7. Computation
%----------------------------------------------------------------

initval;
  l = l_ss;
  k = k_ss;
  c = c_ss;
  w = w_ss;
  r = r_ss;
  u = v_ss;
  v = v_ss;
  rf= 1/beta;

  z = 0;
  sigma_zt = sigma;
  e = 0;
  e_sigma = 0;

  ev = ev_ss;

end;
    
shocks;
  var e = 1;
  var e_sigma = 1;
end;

steady;

set_dynare_seed(2)
stoch_simul(PRUNING, PERIODS = 100000, DROP = 1000, IRF = 0, ORDER = 2) c y i k rf r;