clear
close all

%%%%%%%%%%%%%%%%%%%
% ANDOLFATTO 1996 %
%%%%%%%%%%%%%%%%%%%


% Coefficients and Steady State Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calibration (coefficients)
beta = 0.99;
eta=2.0;
phi1=2.08;
phi2=1.37;
zeta=1.271719032;
theta=0.36;
delta=0.025;
rho=0.95;
epsilon=0.007;
chi=1.011465308;
alpha=0.60;
sigma =0.15;
kappa=0.105;
e=0.165;

% Steady state values
c=0.74;
n=0.57;
l=0.33;
k=10;
v=0.095;
mu=0.157657658;
w=3.345702282;
z=0;

% Displaying the steady states values
ss_values = [c n l k v mu w z];
disp('c  n  l  k  v  mu  w  z');
disp(ss_values);

% Linearized Model in Matrix Form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M1=zeros(6,6);
M2=zeros(6,5);
M3_I=zeros(5,5);
M3_L=zeros(5,5);
M4_I=zeros(5,6);
M4_L=zeros(5,6);
M5=zeros(5,5);


% Static equations
%-----------------

% Equation 1
M1(1,1)=1;
M(1,3)=(l*(n+theta)-1)/(1-l);

M2(1,1)=-theta;
M2(1,2)=theta,
M2(1,3)=1;

% Equation 2
M1(2,1)=-1;
M1(2,6)=1-alpha;

M2(2,2)=-(1-alpha)*n/(1-n);
M2(2,4)=1;

% Equation 3
M1(3,1)=-alpha*(mu*((1-sigma)-1/beta)+phi1*(1-l)^(-eta));
M1(3,3)=(1-alpha)*phi1*(1-l)^(-eta)*c*theta+aloha*phi1*(1-l)^(-eta)*c+alpha*(mu*((1-sigma)-1/beta)+phi1*(1-l)^(-eta));
M1(3,5)=w;
M1(3,6)=-alpha^2*(1-alpha)*sigma*n/(1-n)*mu*c/l;

M2(3,1)=(1-alpha)*phi1*(1-l)^(-eta)*c*theta;
M2(3,2)=-(1-alpha)*phi1*(1-l)^(-eta)*c*theta;
M2(3,4)=alpha*c/l*(1-alpha)*sigma*n*mu/(1-n);

% Equation 4
M1(4,2)=1
M1(4,3)=-1+theta;

M2(4,1)=theta;
M2(4,2)=1-theta;
M2(4,3)=1;

% Equation 5
M1(5,1)=c;
M1(5,3)=-zeta*k^(theta)*(n*l)^(1-theta)*(1-theta);
M1(5,4)=I;
M1(5,6)=kappa*v;

M2(5,1)=zeta*k^(theta)*(n*l)^(1-theta)*theta;
M2(5,2)=zeta*k^(theta)*(n*l)^(1-theta)*(1-theta);
M2(5,3)=zeta*k^(theta)*(n*l)^(1-theta)*;

% Equation 6
M1(6,1)=-1/c;
M2(6,5)=lambda;


% Dynamic Equations
%%%%%%%%%%%%%%%%%%

% Equation d'Euler (P4) - Equation 7
M3_I(1,1)=-(1-beta*(1-delta))*(1-theta);
M3_I(1,2)=(1-beta*(1-delta))*(1-theta);
M3_I(1,3)=(1-beta*(1-delta));
M3_I(1,5)=c*lambda;

M3_L(1,5)=-c*lambda;

M4_I(1,3)=-(1-beta*(1-delta));

M5(1,1)=c*lambda;
M5(1,2)=(1-beta*(1-delta))*(1-theta);
M5(1,3)=(1-beta*(1-delta));

% Equation 8
M3_I(2,1)=-beta*phi1*(1-l)^(-eta)*l*theta;
M3_I(2,2)=beta*phi1*(1-l)^(-eta)*l*theta+beta*mu*(1-alpha)*sigma*n^2*alpha/((1-n)^2);
M3_I(2,3)=-beta*phi1*(1-l)^(-eta)*l;
M3_I(2,4)=-(1-sigma-(1-alpha))*sigma*n/(1-n);
M3_I(2,5)=-beta*phi1*(1-l)^(-eta)*l*c*lambda;

M3_L(2,4)=mu;

M4_I(2,3)=-beta*phi1*(1-l)^(-eta)*l*theta;
M4_I(2,6)=-beta*mu*(1-alpha)*sigma*n*alpha/(1-n);

M5(2,1)=-beta*phi1*(1-l)^(-eta)*l*c*lambda;
M5(2,2)=beta*phi1*(1-l)^(-eta)*l*theta;
M5(2,3)=-beta*phi1*(1-l)^(-eta)*l;
M5(2,4)=-(1-sigma-(1-alpha))*sigma*n/(1-n);
M5(2,5)=beta*mu*(1-alpha)*sigma*n*alpha/(1-n);

% Equation 9
M3_I(3,2)=1;
M3_L(3;2)=-(1-n-sigma*(1-alpha*n)/(1-n));
M4_I(3,6)=sigm*alpha;

% Equation 10
M3_I(4,1)=k;
M3_L(4,1)=-(1-delta);
M4_L(4,4)=I;

% Equation 11
M3_I(5,3)=1;
M3_L(5,3)=-rho;
M5(5,3)=1;


