var c, n, l, k, v, mu, w, z, junk;
varexo epsilon_tilde;
parameters beta, eta, phi1, phi2, zeta, theta, delta, rho, epsilon, chi, alpha, sigma, kappa, e;

beta = 0.99;
eta=2.0;
phi1=2.08;
phi2=1.37;
zeta=0.796562824;
theta=0.36;
delta=0.025;
rho=0.95;
epsilon=0.007;
chi=1.011465308;
alpha=0.60;
sigma =0.15;
kappa=0.105;
e=0.165;

model;
1/c=beta*(theta*exp(z(+1))*zeta*(n(+1)*l(+1)/k(+1))^(1-theta)+1-delta)/c(+1);
phi1*(1-l)^(-eta)=(1-theta)/c*exp(z)*zeta*(k/(n*l))^(theta);
kappa*v/c=mu*alpha*chi*v^(alpha)*((1-n)*e)^(1-alpha);
mu=beta*(phi1*((1-l(+1))^(1-eta)/(1-eta))-phi2/(1-eta)*(1-e)^(1-eta)+1/c(+1)*(1-theta)*exp(z)*zeta*k(+1)^(theta)*(n(+1)*l(+1))^(-theta)*l(+1)+mu(+1)*(1-sigma-(1-alpha)/(1-n(+1))*chi*v(+1)^(alpha)+((1-n(+1))*e)^(1-alpha)));
c+k(+1)+kappa*v=exp(z)*zeta*k^(theta)*(n*l)^(1-theta)+(1-delta)*k;
n(+1)=(1-sigma)*n+chi*v^(alpha)*((1-n)*e)^(1-alpha);
w=(1-alpha)*(1-theta)*exp(z)*zeta*(k/(n*l))^(theta)+alpha*((phi2/(1-eta)*(1-e)^(1-eta)-phi1/(1-eta)*(1-l)^(1-eta))+(1-alpha)*mu/(1-n)*chi*v^(alpha)*((1-n)*e)^(1-alpha))*c/l;
z(+1)=rho*z+epsilon_tilde*epsilon;
junk=0.9*junk(-1);
end;

initval;
c=0.74;
n=0.57;
l=0.33;
k=10;
v=0.095;
mu=0.157657658;
w=3.345702282;
z=0;
junk=0;
end;


shocks;
var epsilon_tilde; stderr 1;
end;


stoch_simul(periods=2100);
                                               
                                                               
          
                                               
                                                        
                                                                        
           
             
                       
