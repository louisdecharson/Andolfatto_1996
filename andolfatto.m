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
sigma=0.15;
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
y=1;
I=1.01-c;
lambda=1/c;


% Displaying the steady states values
ss_values = [c n l k v mu w z];
disp('y c  I n  l  k  v  mu  w  z');
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
M1(1,3)=n*l/(1-l)+theta;

M2(1,1)=theta;
M2(1,2)=-theta,
M2(1,3)=1;

% Equation 2
M1(2,1)=-1;
M1(2,6)=1-alpha;

M2(2,2)=-(1-alpha)*n/(1-n);
M2(2,4)=1;

% Equation 3
M1(3,1)=-alpha*(mu*((1-sigma)-1/beta)+phi1*(1-l)^(-eta));
M1(3,3)=(1-alpha)*phi1*(1-l)^(-eta)*c*theta-alpha*phi1*(1-l)^(-eta)*c+alpha*(mu*((1-sigma)-1/beta)+phi1*(1-l)^(-eta));
M1(3,5)=w;
M1(3,6)=-alpha^2*(1-alpha)*sigma*n/(1-n)*mu*c/l;

M2(3,1)=(1-alpha)*phi1*(1-l)^(-eta)*c*theta;
M2(3,2)=-(1-alpha)*phi1*(1-l)^(-eta)*c*theta;
M2(3,3)=(1-alpha)*phi1*(1-l)^(-eta)*c;
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
M2(5,3)=zeta*k^(theta)*(n*l)^(1-theta);

% Equation 6
M1(6,1)=-1;
M2(6,5)=1;


% Dynamic Equations
%%%%%%%%%%%%%%%%%%

% Equation d'Euler (P4) - Equation 7
M3_I(1,1)=-(1-beta*(1-delta))*(1-theta);
M3_I(1,2)=(1-beta*(1-delta))*(1-theta);
M3_I(1,3)=(1-beta*(1-delta));
M3_I(1,5)=c*lambda;

M3_L(1,5)=-c*lambda;

M4_I(1,3)=-(1-beta*(1-delta))*(1-theta);

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

M3_L(3,2)=-(1-n-sigma*(1-alpha*n))/(1-n);

M4_I(3,6)=sigma*alpha;

% Equation 10
M3_I(4,1)=k;
M3_L(4,1)=-(1-delta)*k;
M4_L(4,4)=I;

% Equation 11
M3_I(5,3)=1;
M3_L(5,3)=-rho;
M5(5,3)=1;


% Reduced Form
%%%%%%%%%%%%%%

M1_inv=inv(M1);

L0=M3_I-M4_I*M1_inv*M2;
L1=M4_L*M1_inv*M2-M3_L;
L2=M5;

W=inv(L0)*L1;
Q=inv(L0)*L2;

% Eigenvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[PP,DD]=eig(W);

% PP is the matrix of the eigen vectors
% DD is the diagonal matrix of the eigen values so that W*PP=PP*DD

[llambda,kk]=sort(diag(DD));  %kk keeps the ranking of the eigenvalues 

P1=PP(:,kk);
P1I=inv(P1);
MU1=P1I*W*P1;

[abf abf] = size(MU1);

ab=0;
for i = 1:abf;
if abs(MU1(i,i)) < 1,
ab=ab+1;
else;
end;
end;

af = abf- ab;

for i = 1:abf;
if abs(MU1(i,i)) == 1,
disp('Unit root')
else;
end;
end;

disp('backward    ')
disp(ab)
disp('forward')
disp(af)

disp('Eigenvalues')
disp(' ')
disp(diag(MU1));

% Saddle Path Condition 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coordonnées / variables backward
P1IB=[P1I(4:5,1:3)];

% coordonnées / variables forward
P1IF=[-P1I(4:5,4:5)];

% condition initiale
GG=inv(P1IF)*P1IB;


% Policy rules 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dynamique variable backward
WB=[W(1:3,1:3)];
WF=[W(1:3,4:5)];
PIB=WB+(WF*GG);

%dynamique variable de contrôle
M2B=[M2(1:6,1:3)];
M2F=[M2(1:6,4:5)];
PIC=M1\(M2B+M2F*GG);

% Regle de decisions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('Policy rules(k,n,z) state variables') 
disp(' ')
disp(PIB)


disp('Policy rules (c,y,l,I,w,v) control variables')
disp(' ')
disp(PIC)

% Impulse Response Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrep=100;

CHOC=[0 ; 0 ; 1]

for j=1:nrep;
RC=(PIB^(j-1))*CHOC;
RCK(j)=RC(1);
RCN(j)=RC(2);
RCZ(j)=RC(3);
end;

for j=1:nrep;
RC=(PIC*PIB^(j-1))*CHOC;
 RCC(j)=RC(1);
 RCY(j)=RC(2);
 RCL(j)=RC(3);
 RCI(j)=RC(4);
 RCW(j)=RC(5);
 RCV(j)=RC(6);
end;

figure
subplot(221),plot(RCK(1:100))
title('Capital stock (K)')
xlabel('quarters')
ylabel('% Dev.   ')

subplot(222),plot(RCZ(1:100))
title('Productivity (Z)')
xlabel('quarters')
ylabel('% Dev.   ')

subplot(223),plot(RCW(1:100))
title('Wages (W)')
xlabel('quarters')
ylabel('% Dev.   ')

subplot(224),plot(RCN(1:100))
title(' Employment rate (N)')
xlabel('quarters')
ylabel('% Dev.   ')


figure
subplot(221),plot(RCC(1:100))
title('Consumption')
xlabel('quarters')
ylabel('% Dev.   ')

subplot(222),plot(RCY(1:100))
title('Output')
xlabel('quarters')
ylabel('% Dev.   ')


subplot(223),plot(RCL(1:100))
title('Hours')
xlabel('quarters')
ylabel('% Dev.   ')

subplot(224),plot(RCI(1:100))
title(' Investment')
xlabel('quarters')
ylabel('% Dev.   ')


figure
subplot(221),plot(RCV(1:100))
title(' New Jobs')
xlabel('quarters')
ylabel('% Dev.   ')


% Stochastic simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nsimul=1000;

nlong=160;


for j=1:nsimul;

disp('simulation')
disp(j)

% simulation des jiemes parties cycliques 

aleaa(1:nlong,j)=randn(nlong,1);

% tirage des innovations

for i=1:nlong;
epsa(i)= aleaa(i,j) * epsilon;
end;


% construction des chocs CHT

CHT(1)=epsa(1);

for i=2:nlong;
CHT(i)=rho*CHT(i-1)+epsa(i);
end;

% initialisation de la partie cyclique du capital et de l'emploi

KNC=zeros(2,nlong);
KNC(1,1)=0;
KNC(2,1)=0;

% parties cycliques

for i=2:nlong;
KNC(:,i)=PIB(1:2,1:2)*KNC(:,i-1)+PIB(1,3)*CHT(i-1);

end;

for i=1:nlong;

CC(i)=PIC(1,1:2)*KNC(:,i)+PIC(1,3)*CHT(i);

YC(i)=PIC(2,1:2)*KNC(:,i)+PIC(2,3)*CHT(i);

LC(i)=PIC(3,1:2)*KNC(:,i)+PIC(3,3)*CHT(i);

IC(i)=PIC(4,1:2)*KNC(:,i)+PIC(4,3)*CHT(i);

WC(i)=PIC(5,1:2)*KNC(:,i)+PIC(5,3)*CHT(i);

VC(i)=PIC(6,1:2)*KNC(:,i)+PIC(6,3)*CHT(i);

end;


% Business Cycles Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[CHPT,CHP]=hpfilter(CC,1600);
[YHPT,YHP]=hpfilter(YC,1600);
[LHPT,LHP]=hpfilter(LC,1600);
[IHPT,IHP]=hpfilter(IC,1600);
[WHPT,WHP]=hpfilter(WC,1600);
[VHPT,VHP]=hpfilter(VC,1600);



ETCHP(j)=std(CHP(1:nlong));
ETYHP(j)=std(YHP(1:nlong));
ETLHP(j)=std(LHP(1:nlong));
ETIHP(j)=std(IHP(1:nlong));
ETWHP(j)=std(WHP(1:nlong));
ETVHP(j)=std(VHP(1:nlong));


RHO=corrcoef(YHP,CHP);
RHOCHP(j)=RHO(1,2);
RHO=corrcoef(YHP,LHP);
RHOLHP(j)=RHO(1,2);
RHO=corrcoef(YHP,IHP);
RHOIHP(j)=RHO(1,2);
RHO=corrcoef(YHP,WHP);
RHOWHP(j)=RHO(1,2);
RHO=corrcoef(YHP,VHP);
RHOVHP(j)=RHO(1,2);


xxC=CHP(2:nlong);
yyC=CHP(1:nlong-1);
RHO=corrcoef(xxC,yyC);
RHOCCHP(j)=RHO(1,2);

xxY=YHP(2:nlong);
yyY=YHP(1:nlong-1);
RHO=corrcoef(xxY,yyY);
RHOYYHP(j)=RHO(1,2);

xxN=LHP(2:nlong);
yyN=LHP(1:nlong-1);
RHO=corrcoef(xxN,yyN);
RHOLLHP(j)=RHO(1,2);


xxI=IHP(2:nlong);
yyI=IHP(1:nlong-1);
RHO=corrcoef(xxI,yyI);
RHOIIHP(j)=RHO(1,2);


xxP=WHP(2:nlong);
yyP=WHP(1:nlong-1);
RHO=corrcoef(xxP,yyP);
RHOWWHP(j)=RHO(1,2);

xxP=VHP(2:nlong);
yyP=VHP(1:nlong-1);
RHO=corrcoef(xxP,yyP);
RHOVVHP(j)=RHO(1,2);

end;


mETCHP=mean(ETCHP);
mETYHP=mean(ETYHP);
mETLHP=mean(ETLHP);
mETIHP=mean(ETIHP);
mETWHP=mean(ETWHP);
mETVHP=mean(ETVHP);



mRHOCHP=mean(RHOCHP);
mRHOLHP=mean(RHOLHP);
mRHOIHP=mean(RHOIHP);
mRHOWHP=mean(RHOWHP);
mRHOVHP=mean(RHOVHP);


mRHOCCHP=mean(RHOCCHP);
mRHOYYHP=mean(RHOYYHP);
mRHOLLHP=mean(RHOLLHP);
mRHOIIHP=mean(RHOIIHP);
mRHOWWHP=mean(RHOWWHP);
mRHOVVHP=mean(RHOVVHP);

mET1=[mETCHP mETYHP mETLHP mETIHP mETWHP mETVHP]./mETYHP;
mRHO1=[mRHOCHP  1  mRHOLHP mRHOIHP mRHOWHP mRHOVHP];
mARHO1=[mRHOCCHP mRHOYYHP mRHOLLHP  mRHOIIHP mRHOWWHP mRHOVVHP];

disp('variable order: C - Y - L - I - W - V')
disp(' ')

disp('relative standard deviation')
disp(' ')
disp(mET1)


disp('correlation')
disp(' ')
disp(mRHO1)


disp('serial correlation')
disp(' ')
disp(mARHO1)
