clear
close all

%%%%%%%%%%%%%%%%%%%
% ANDOLFATTO 1996 %
%%%%%%%%%%%%%%%%%%%


% Coefficients and Steady State Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calibration (coefficients)

% We know :
beta=0.99;
delta=0.025;
theta=0.36;

% zeta is chosen to normalize steady state output ->
% zeta=k^(-theta)*(n*l)^(theta-1)
% Then k is computed from (P1') (p.120)
y=1;
z=0;
k=theta/(1/beta-1+delta);

% We know n and l
n=0.57;
l=0.33;
zeta=k^(-theta)*(n*l)^(theta-1);

% We know that psi=kappa*v=0.01, we compute c from P5'
psi=0.01;
c=1-psi-delta*k

% We compute kappa and v from equation 28 (page 121)
sigma=0.15;
q=0.90;
v=sigma*n/q;
kappa=psi/v;

% We compute chi and mu from P3' and from the matching function
% definition
% mu = psi/(c*sigma*alpha*n)
e=l/2;
alpha=0.6;
chi=sigma*n/(v^(alpha)*((1-n)*e)^(1-alpha));
mu=psi/(c*sigma*alpha*n);



% page 121
eta=2;
rho=0.95;
epsilon=0.007;

% For wages, we use equation 27
w=((1-theta)-(1-(1-sigma)*beta)*psi/(beta*sigma))/(n*l);


% Other coefficients :
phi1=2.08;


I=delta*k;
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
M2(1,2)=-theta;
M2(1,3)=1;

% Equation 2
M1(2,1)=-1;
M1(2,6)=1-alpha;

M2(2,2)=-(1-alpha)*n/(1-n);
M2(2,4)=1;

% Equation 3
M1(3,1)=-alpha*c/l*(mu*((1-sigma)-1/beta)+phi1*l*(1-l)^(-eta));
M1(3,3)=(1-alpha)*phi1*(1-l)^(-eta)*c*theta-alpha*phi1*(1-l)^(-eta)*c+alpha*c/l*(mu*((1-sigma)-1/beta)+phi1*l*(1-l)^(-eta));
M1(3,5)=w;
M1(3,6)=-alpha^2*(1-alpha)*sigma*n/(1-n)*mu*c/l;

M2(3,1)=(1-alpha)*phi1*(1-l)^(-eta)*c*theta;
M2(3,2)=-(1-alpha)*phi1*(1-l)^(-eta)*c*theta+mu*c/l*(alpha^2*(1-alpha)*sigma*n^2)/((1-n)^2);
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

nsimul=100;

nlong=160;


for j=1:nsimul;

    disp('simulation')
    disp(j)

    % simulation des jiemes parties cycliques 

    for i=1:nlong
        aleaa(i,j)=binornd(1,0.5);
        if aleaa(i,j)==0;
            aleaa(i,j)=-1;
        end;
    end;


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
        KNC(1,i)=PIB(1,1:2)*KNC(:,i-1)+PIB(1,3)*CHT(i-1);
        KNC(2,i)=PIB(2,1:2)*KNC(:,i-1)+PIB(2,3)*CHT(i-1);
    end;

    for i=1:nlong;

        CC(i)=PIC(1,1:2)*KNC(:,i)+PIC(1,3)*CHT(i);

        YC(i)=PIC(2,1:2)*KNC(:,i)+PIC(2,3)*CHT(i);

        LC(i)=PIC(3,1:2)*KNC(:,i)+PIC(3,3)*CHT(i);

        IC(i)=PIC(4,1:2)*KNC(:,i)+PIC(4,3)*CHT(i);

        WC(i)=PIC(5,1:2)*KNC(:,i)+PIC(5,3)*CHT(i);

        VC(i)=PIC(6,1:2)*KNC(:,i)+PIC(6,3)*CHT(i);

        % Total hours
        LNC(i)=KNC(2,i)+LC(i);
        % Wage Bill 
        WLNC(i)=LNC(i)+WC(i);
        
        % Labor share
        LSC(i)=WLNC(i)-YC(i);
        
        % Productivity
        PC(i)=YC(i)-LNC(i);
        
        % Unemployment
        UC(i)=-n/(1-n)*KNC(2,i);
    end;

 
        
    % Business Cycles Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [CHPT,CHP]=hpfilter(CC,1600);
    [YHPT,YHP]=hpfilter(YC,1600);
    [LHPT,LHP]=hpfilter(LC,1600);
    [IHPT,IHP]=hpfilter(IC,1600);
    [WHPT,WHP]=hpfilter(WC,1600);
    [VHPT,VHP]=hpfilter(VC,1600);

    %[KHPT,KHP]=hpfilter(KNC(1,:),1600);
    [NHPT,NHP]=hpfilter(KNC(2,:),1600);

    [LNHPT,LNHP]=hpfilter(LNC,1600);
    [WLNHPT,WLNHP]=hpfilter(WLNC,1600);
    [LSHPT,LSHP]=hpfilter(LSC,1600);
    [PHPT,PHP]=hpfilter(PC,1600);
    [UHPT,UHP]=hpfilter(UC,1600);

    % Standard Deviation

    ETCHP(j)=std(CHP(1:nlong));
    ETYHP(j)=std(YHP(1:nlong));
    ETLHP(j)=std(LHP(1:nlong));
    ETIHP(j)=std(IHP(1:nlong));
    ETWHP(j)=std(WHP(1:nlong));
    ETVHP(j)=std(VHP(1:nlong));
    
    %ETKHP(j)=std(KHP(1:nlong));
    ETNHP(j)=std(NHP(1:nlong));

    ETLNHP(j)=std(LNHP(1:nlong));
    ETWLNHP(j)=std(WLNHP(1:nlong));
    ETLSHP(j)=std(LSHP(1:nlong));
    ETPHP(j)=std(PHP(1:nlong));
    ETUHP(j)=std(UHP(1:nlong));
    
    % Correlation

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

    %RHO=corrcoef(YHP,KHP);
    %RHOKHP(j)=RHO(1,2);
    RHO=corrcoef(YHP,NHP);
    RHONHP(j)=RHO(1,2);

    RHO=corrcoef(YHP,LNHP);
    RHOLNHP(j)=RHO(1,2);

    RHO=corrcoef(YHP,WLNHP);
    RHOWLNHP(j)=RHO(1,2);

    RHO=corrcoef(YHP,LSHP);
    RHOLSHP(j)=RHO(1,2);

    RHO=corrcoef(YHP,PHP);
    RHOPHP(j)=RHO(1,2);
    
    RHO=corrcoef(YHP,UHP);
    RHOUHP(j)=RHO(1,2);
    
    %% Cross Correlation
    
    % of Total Hours with Productivity and Real Wage
    % And then of Unemployment with Unemployment and Vacancies
    
    for kk=1:4
        
        % Total Hours with Productivity
        xxLN=LNHP(1+kk:nlong);
        yyP=PHP(1:nlong-kk);
        RHO=corrcoef(xxLN,yyP);
        CrossRHO_P(j,5-kk)=RHO(1,2);
        
        RHO=corrcoef(LNHP,PHP);
        CrossRHO_P(j,5)=RHO(1,2);
        
        xxLN=LNHP(1:nlong-kk);
        yyP=PHP(1+kk:nlong);
        RHO=corrcoef(xxLN,yyP);
        CrossRHO_P(j,5+kk)=RHO(1,2);
        
        % Total Hours with Real Wage
        xxLN=LNHP(1+kk:nlong);
        yyW=WHP(1:nlong-kk);
        RHO=corrcoef(xxLN,yyW);
        CrossRHO_W(j,5-kk)=RHO(1,2);
        
        RHO=corrcoef(LNHP,WHP);
        CrossRHO_W(j,5)=RHO(1,2);
        
        xxLN=LNHP(1:nlong-kk);
        yyW=WHP(1+kk:nlong);
        RHO=corrcoef(xxLN,yyW);
        CrossRHO_W(j,5+kk)=RHO(1,2);
        
        % Unemployment with Unemployment
        xxU=UHP(1+kk:nlong);
        yyU=UHP(1:nlong-kk);
        RHO=corrcoef(xxU,yyU);
        CrossRHO_U(j,5-kk)=RHO(1,2);
        
        RHO=corrcoef(UHP,UHP);
        CrossRHO_U(j,5)=RHO(1,2);
        
        xxU=UHP(1:nlong-kk);
        yyU=UHP(1+kk:nlong);
        RHO=corrcoef(xxU,yyU);
        CrossRHO_U(j,5+kk)=RHO(1,2);
        
        % Unemployment with Vacancies
        xxU=UHP(1+kk:nlong);
        yyV=VHP(1:nlong-kk);
        RHO=corrcoef(xxU,yyV);
        CrossRHO_V(j,5-kk)=RHO(1,2);
        
        RHO=corrcoef(UHP,VHP);
        CrossRHO_V(j,5)=RHO(1,2);
        
        xxU=UHP(1:nlong-kk);
        yyV=VHP(1+kk:nlong);
        RHO=corrcoef(xxU,yyV);
        CrossRHO_V(j,5+kk)=RHO(1,2);
        
    end
    
    
    % Serial Correlation

    xxC=CHP(2:nlong);
    yyC=CHP(1:nlong-1);
    RHO=corrcoef(xxC,yyC);
    RHOCCHP(j)=RHO(1,2);

    xxY=YHP(2:nlong);
    yyY=YHP(1:nlong-1);
    RHO=corrcoef(xxY,yyY);
    RHOYYHP(j)=RHO(1,2);

    xxL=LHP(2:nlong);
    yyL=LHP(1:nlong-1);
    RHO=corrcoef(xxL,yyL);
    RHOLLHP(j)=RHO(1,2);


    xxI=IHP(2:nlong);
    yyI=IHP(1:nlong-1);
    RHO=corrcoef(xxI,yyI);
    RHOIIHP(j)=RHO(1,2);


    xxW=WHP(2:nlong);
    yyW=WHP(1:nlong-1);
    RHO=corrcoef(xxW,yyW);
    RHOWWHP(j)=RHO(1,2);

    xxV=VHP(2:nlong);
    yyV=VHP(1:nlong-1);
    RHO=corrcoef(xxV,yyV);
    RHOVVHP(j)=RHO(1,2);

    %xxK=KHP(2:nlong);
    %yyK=KHP(1:nlong-1);
    %RHO=corrcoef(xxK,yyK);
    %RHOKKHP(j)=RHO(1,2);

    xxN=NHP(2:nlong);
    yyN=NHP(1:nlong-1);
    RHO=corrcoef(xxN,yyN);
    RHONNHP(j)=RHO(1,2);

    xxLN=LNHP(2:nlong);
    yyLN=LNHP(1:nlong-1);
    RHO=corrcoef(xxLN,yyLN);
    RHOLLNNHP(j)=RHO(1,2);

    xxWLN=WLNHP(2:nlong);
    yyWLN=WLNHP(1:nlong-1);
    RHO=corrcoef(xxWLN,yyWLN);
    RHOWWLLNNHP(j)=RHO(1,2);

    xxLS=LSHP(2:nlong);
    yyLS=LSHP(1:nlong-1);
    RHO=corrcoef(xxLS,yyLS);
    RHOLLSSHP(j)=RHO(1,2);
    
    xxP=PHP(2:nlong);
    yyP=PHP(1:nlong-1);
    RHO=corrcoef(xxP,yyP);
    RHOPPHP(j)=RHO(1,2);
    
    xxU=UHP(2:nlong);
    yyU=UHP(1:nlong-1);
    RHO=corrcoef(xxU,yyU);
    RHOUUHP(j)=RHO(1,2);

end;

% We take the mean over the simulations

% Relative Standard Deviation
mETCHP=mean(ETCHP);
mETYHP=mean(ETYHP);
mETLHP=mean(ETLHP);
mETIHP=mean(ETIHP);
mETWHP=mean(ETWHP);
mETVHP=mean(ETVHP);
%mETKHP=mean(ETKHP);
mETNHP=mean(ETNHP);
mETLNHP=mean(ETLNHP);
mETWLNHP=mean(ETWLNHP);
mETLSHP=mean(ETLSHP);
mETPHP=mean(ETPHP);
mETUHP=mean(ETUHP);

% Correlation 
mRHOCHP=mean(RHOCHP);
mRHOLHP=mean(RHOLHP);
mRHOIHP=mean(RHOIHP);
mRHOWHP=mean(RHOWHP);
mRHOVHP=mean(RHOVHP);
%mRHOKHP=mean(RHOKHP);
mRHONHP=mean(RHONHP);
mRHOLNHP=mean(RHOLNHP);
mRHOWLNHP=mean(RHOWLNHP);
mRHOLSHP=mean(RHOLSHP);
mRHOPHP=mean(RHOPHP);
mRHOUHP=mean(RHOUHP);

% Cross Correlation
mCrossRHO_W=mean(CrossRHO_W,1);
mCrossRHO_P=mean(CrossRHO_P,1);
mCrossRHO_U=mean(CrossRHO_U,1);
mCrossRHO_V=mean(CrossRHO_V,1);



% serial CORRELATION 
mRHOCCHP=mean(RHOCCHP);
mRHOYYHP=mean(RHOYYHP);
mRHOLLHP=mean(RHOLLHP);
mRHOIIHP=mean(RHOIIHP);
mRHOWWHP=mean(RHOWWHP);
mRHOVVHP=mean(RHOVVHP);

%mRHOKKHP=mean(RHOKKHP);
mRHONNHP=mean(RHONNHP);

mRHOLLNNHP=mean(RHOLLNNHP);
mRHOWWLLNNHP=mean(RHOWWLLNNHP);
mRHOLLSSHP=mean(RHOLLSSHP);
mRHOPPHP=mean(RHOPPHP);
mRHOUUHP=mean(RHOUUHP);


% mET1=[mETCHP mETYHP mETLHP mETIHP mETWHP mETVHP mETNHP]./mETYHP;
% %mETLNHP mETWLNHP mETRWHP
% mRHO1=[mRHOCHP  1  mRHOLHP mRHOIHP mRHOWHP mRHOVHP mRHONHP mRHOLNHP mRHOWLNHP mRHORWHP];
% mARHO1=[mRHOCCHP mRHOYYHP mRHOLLHP  mRHOIIHP mRHOWWHP mRHOVVHP mRHONNHP mRHOLLNNHP mRHOWWLLNNHP mRHORRWWHP];


% disp('variable order: C - Y - L - I - W - V  - N - LN - WLN - RW')
% disp(' ')

% disp('relative standard deviation')
% disp(' ')
% disp(mET1)


% disp('correlation')
% disp(' ')
% disp(mRHO1)


% disp('serial correlation')
% disp(' ')
% disp(mARHO1)


disp('Consumption (C)')
disp( 'rel St.Dev   Corr      ser corr')
disp([mETCHP/mETYHP mRHOCHP mRHOCCHP])

disp('Investment (I)')
disp( 'rel St.Dev   Corr      ser corr')
disp([mETIHP/mETYHP mRHOIHP mRHOIIHP])

disp(' ')
disp(' ')
disp('Total Hours (L*N)')
disp( 'rel St.Dev   Corr      ser corr')
disp([mETLNHP/mETYHP mRHOLNHP mRHOLLNNHP])

disp('Employment (N)')
disp( 'rel St.Dev   Corr      ser corr')
disp([mETNHP/mETYHP mRHOIHP mRHOIIHP])

disp('Hours/worker (L)')
disp( 'rel St.Dev   Corr      ser corr')
disp([mETLHP/mETYHP mRHOLHP mRHOLLHP])

disp(' ')
disp(' ')
disp('Wage bill (WLN)')
disp( 'rel St.Dev   Corr      ser corr')
disp([mETWLNHP/mETYHP mRHOWLNHP mRHOWWLLNNHP])


disp('Labor Share (LS)')
disp( 'rel St.Dev   Corr      ser corr')
disp([mETLSHP/mETYHP mRHOLSHP mRHOLLSSHP])


disp(' ')
disp(' ')
disp('Productivity (P)')
disp( 'rel St.Dev   Corr      ser corr')
disp([mETPHP/mETYHP mRHOPHP mRHOPPHP])

disp('Real Wage (W)')
disp( 'rel St.Dev   Corr      ser corr')
disp([mETWHP/mETYHP mRHOWHP mRHOWWHP])

disp(' ')
disp(' ')
disp('Unemployment (U)')
disp( 'rel St.Dev   Corr      ser corr')
disp([mETUHP/mETYHP mRHOUHP mRHOUUHP])
disp('Vacancies (U)')
disp( 'rel St.Dev   Corr      ser corr')
disp([mETVHP/mETYHP mRHOVHP mRHOVVHP])


disp(' ')
disp(' ')
disp(' ')
disp('Cross Correlations of Hours with Productivity and the Real Wage')
disp('x(t-4)  x(t-3)  x(t-2)  x(t-1)  x(t)  x(t+1)  x(t+2)  x(t+3)  x(t+4)')
disp(mCrossRHO_P)
disp(mCrossRHO_W)

disp('Cross Correlations of Unemployment with Unemployment and the Vacancies')
disp('x(t-4)  x(t-3)  x(t-2)  x(t-1)  x(t)  x(t+1)  x(t+2)  x(t+3)  x(t+4)')
disp(mCrossRHO_U)
disp(mCrossRHO_V)

