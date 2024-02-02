%% Formulated and Encoded by Aninda Pal
% Organization :: IIT Kharagpur, Departmentof Ocean Engg & Naval Architecture
% cite using DOI :: 10.1098/rspa.2023.0707
%% Problem initialization
% Functions required : makefunTeich, Teicheq & build_F
% comparing Teich low intensity
clc
clear all
syms phi psi xia xis rho c m omega Cd a1 gamma alpha t_d x p
syms psi_1(t) psi_2(t) P_ior P_io P_o

% Nomencleture
% P_io= incident overpressure; P_ior= incident overpressure ratio;
% P_o= atmospheric pressure
tic
p_x=P_io; % incident overpressure
% P_o=101325; % atmospheric pressure
P_ior= p_x/P_o; % incident over pressure ratio
psi=P_ior*psi_1(t);

% Reflection coefficient :: For comparison
C_r=((3*gamma-1)*psi+4*gamma)/((1*gamma-1)*psi+2*gamma); % coefficient of reflection
phi=P_ior*psi_2(t)/(1+P_ior*psi_1(t));

% psi
psid= diff(psi,t);
psidd= diff(psid,t);
% phi
phid= diff(phi,t);
phidd= diff(phid,t);

% Simplification stage
A= sqrt(1+(gamma+1)/(2*gamma)*(psi))*a1;
B= 2*psi/(4*gamma+2*gamma*psi);
C= phi/psi*sqrt((2*gamma+(gamma-1)*psi)/(2*gamma+(gamma+1)*phi)*(1+psi));
Ex_dot= A*B*(-C+1);% denotes the velocity

% Derivatives
Ex_ddot= diff(Ex_dot,t);
Ex_dddot= diff(Ex_ddot,t);
psi_1dot=diff(psi_1(t),t);
psi_2dot=diff(psi_2(t),t);

% Backpressure term
ad=rho*((gamma+1)/4*Ex_dot+sqrt(((gamma+1)/4)^2*Ex_dot^2+a1^2))*Ex_dot; % backpressure damping

% Resistive part of Governing differential equation
EL= m*(Ex_dddot)+2*m*omega*xis*(Ex_ddot)+1*diff(ad,t)+m*omega^2*Ex_dot;% For nonlinear damping
% EL= m*(Ex_dddot)+2*m*omega*xis*(Ex_ddot)+Cd*rho*a1*(Ex_ddot)+m*omega^2*Ex_dot; % For Linear damping
% EL= m*(Ex_dddot)+2*m*omega*xis*(Ex_ddot)+0*(Ex_ddot)+m*omega^2*Ex_dot; % For Constant damping

% Exciting part of Governing differential equation
ER=P_o*P_ior*(psi_1dot+psi_2dot);


% Parameters :: we change accordingly
omg=200;% frequency in hz
mass=50;% areal mass
cod=1; % drag coefficient
sd=0.03;% structural damping
inop=6829e3; % incident over pressure
timd=2.33e-3;   % time duration
Pz=101325;% atmospheric pressure
ss= 340; % sound speed at stp 
dc= 7.99; % decay coefficient 
gam=1.4;% specific heat ratio
ad=1.225;% air density 

% Find the initial value problem
syms AA AB
PR=(1-t/t_d)*exp(-alpha*t/t_d);% profile
% PR=exp(-alpha*t/t_d);% profile
Eq=subs(Ex_dot,[psi_1 psi_2 t],[PR AA 0])==0;
Eq1=subs(Eq,[P_o P_io Cd a1 alpha gamma m omega rho t_d xis],[Pz inop cod ss dc gam mass omg*2*pi ad timd sd]);
Eq2=double(solve(Eq1,AA));
Line1= 'Reflected over-pressure profile amplitude is ';
disp([Line1, num2str(Eq2),' units']);
disp(['When the amplitude of incident overpressure profile is 1 units']);
C_rval=subs(C_r,[psi_1 P_o P_io Cd a1 alpha gamma m omega rho t_d xis t],[PR Pz inop cod ss dc gam mass omg*2*pi ad timd sd 0]);
double(C_rval);
Line2= 'Reflection coefficient is ';
disp([Line2, num2str(double(C_rval))]);

ip=PR;
ip=subs(ip,[alpha t_d],[dc timd]);
Eq3=Ex_ddot;
Eq4=subs(Eq3,[psi_1(t) diff(psi_2,t)],[PR AB]);
Eq5=Eq4==P_o*P_ior*(1+Eq2)/mass;
Eq6=solve(Eq5,AB);
Eq7=subs(Eq6,[psi_1 P_o P_io Cd a1 alpha gamma m omega rho t_d xis t],[PR Pz inop cod ss dc gam mass omg*2*pi ad timd sd 0]);
Eq8=double(subs(Eq7,[psi_2(0)],[Eq2]));
syms x(t) p1 p2 p3 
EQ= (EL-ER)==0;
N= EL-ER;
% p1 is acceleration p2 is velocity and p3 is dispalacement
EQ2= subs(N,[diff(diff(psi_2,t)) diff(psi_2,t) psi_2],[p1 p2 p3]);
EQ3=solve(EQ2,p1);
EQ4=subs(EQ3,[psi_1 P_o P_io Cd a1 alpha gamma m omega rho t_d xis],[PR Pz inop cod ss dc gam mass omg*2*pi ad timd sd]);
EQ4=subs(EQ4,[psi_1 P_o P_io Cd a1 alpha gamma m omega rho t_d xis],[PR Pz inop cod ss dc gam mass omg*2*pi ad timd sd]);
EQ5= matlabFunction(EQ4);
[tg1, yg1] = ode23(@(t,y)PG_eq(t,y,EQ5),[0 0.01],[Eq2 Eq8]);
% Plotting
figure(4),clf
pi=double(subs(ip,t,tg1));
hold on
plot(tg1, pi,'k'); % incident profile
plot(tg1, yg1(:,1), 'r-'); % reflected profile
legend({'incident priofile','reflected profile'}, 'Location', 'Best');
title('Wave profiles');
xlabel("t (sec)");
xlim([-0.001 tg1(end)])
grid on
% Refining :: redundant step
si1=subs(PR,[t_d alpha],[timd dc]);
si1=double(subs(si1,t,tg1));
si2=yg1(:,1);
% Find response
syms r1 r2
Eqm=subs(Ex_dot,[P_o P_io a1 gamma],[Pz inop ss gam]);
Eqm1=subs(Eqm,[psi_1 psi_2],[r1 r2]);%si1 si2
Eqm2=double(subs(Eqm1,{r1 r2},{si1 si2})); % gives velocity
dispval=cumtrapz(tg1,Eqm2); % gives displacement

figure(7),clf
plot(tg1,dispval,'m');
xlim([-0.001 1])
title('Displacement response');
xlim([-0.001 tg1(end)])
xlabel("t (sec)");
ylabel("Displacement (m/s)");
grid on