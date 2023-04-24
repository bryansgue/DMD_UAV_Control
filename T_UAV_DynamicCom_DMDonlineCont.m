%******************************************************************************************************************
%************************************ SEGUIMIENTO DE TRAYECTORIA **************************************************
%************************************* ROBOT MANIPULADOR AÉREO *****************************************************
%******************************************************************************************************************
clc; clear all; close all; warning off % Inicializacion
f = 30;
ts = 1/f;       % Tiempo de muestreo
tfin = 120;      % Tiempo de simulación
t = 0:ts:tfin;
a = 0.0; 
b = 0.0;
c = 0.0;
L = [a, b, c];

%% Variables definidas por la TRAYECTORIA y VELOCIDADES deseadas
[xd, yd, zd, psid, xdp, ydp, zdp, psidp] = TrayectoriasDMD(2,t);
%% GENERALIZED DESIRED SIGNALS
hd = [xd; yd; zd; psid];
hd_p = [xdp;ydp;zdp;psidp];

%% a) Posiciones iniciales del UAV
xu(1) = 0; 
yu(1) = 0; 
zu(1) = 5; 
psi(1)= 0;
h=[xu(1);yu(1);zu(1);psi(1)];

%% Matriz A & B del sistema dinamico DMD Offline (Continuo)
load("A_B_values.mat"); 

load("chi_values.mat");
chi_real = chi';

%% Values init for DMD ONLINE (Discreto)
load("G&P_DMDonline_values_init.mat");
Ae = G(:,1:4);
Be = G(:,5:end);

%% Velocidad inicial real del UAV
v = [0;0;0;0];
v_est = v(:,1);
sample = 0;
A_E_P = reshape(Ae,16,1)
B_E_P = reshape(Be,16,1) 

A_c = (Ae-eye(4))/ts;
B_c = Be/ts;


A_c_aux = (Ae-eye(4))/ts;
B_c_aux = Be/ts;



%% Windowed DMD Online Value
m=6;
%******************************************************************************************************************
%***************************************** CONTROLADOR ***********************************************************
%*****************************************************************************************************************
disp('Empieza el programa')
F_extern = zeros(4, length(t));
t_aux = (t>=20 & t<=60) | (t>=80 & t<=100);

F_extern(1,t_aux) = -0.4;
F_extern(2,t_aux) = 0.5;
F_extern(3,t_aux) = -0.1;

%% Filter 
% Filter force design
kp_f = 1;
wn_f = sqrt(kp_f);
kv_f = 1.1*1*wn_f;

Filter = tf([1], [1 kv_f kp_f]);
Filter_d = c2d(Filter, ts);
[num1d_filter, den1d_filter] = tfdata(Filter_d,'v');

% Filter coeficients
A_filter = num1d_filter(2:end);
B_filter = den1d_filter(2:end);

% Filter signals
% Force Filter
ul_memory = zeros(length(den1d_filter(2:end)),1);
um_memory = zeros(length(den1d_filter(2:end)),1);
un_memory = zeros(length(den1d_filter(2:end)),1);
w_memory = zeros(length(den1d_filter(2:end)),1);

ul_filter =  zeros(length(num1d_filter(2:end)),1);
um_filter =  zeros(length(num1d_filter(2:end)),1);
un_filter =  zeros(length(num1d_filter(2:end)),1);
w_filter =  zeros(length(num1d_filter(2:end)),1);

for k=1:length(t)-1
tic
%% 1) LEY DE CONTROL
vc(:,k) = Vc_UAV(hd_p(:,k),hd(:,k),xu(k),yu(k),zu(k),psi(k)); 
ul(k)=vc(1,k); um(k)=vc(2,k); un(k)=vc(3,k); w(k)=vc(4,k);
if k==1
    ulp = ul(k)/ts; ump = um(k)/ts;  unp= un(k)/ts; wp= w(k)/ts;
else
    ulp = (ul(k)-ul(k-1))/ts; ump = (um(k)-um(k-1))/ts;
    unp = (un(k)-un(k-1))/ts; wp = (w(k)-w(k-1))/ts;
end
 vcp = [ulp;ump;unp;wp];
vcp = [0;0;0;0];
%% DYAMIC COMPENSATION
vref(:,k) = dynamicComDMD_online(A_c_aux, B_c_aux, vcp, vc(:,k), v(:,k), 3, 4); %2 6
vref(:,k) = vc(:,k);
%% 2) DINAMICA DEL UAV (VELOCIDAD Y POSICION)
%[v_real(:, k+1),Tu(:,k)] = dyn_model_adapUAV(chi_real, v_real(:,k), vref(:,k), psi(k), L, ts, k, F_extern(:, k));
[v(:, k+1),Tu(:,k)] = DMD_dymamic_system(A,B,v(:,k), vref(:,k),ts,F_extern(:, k));

% Integracion numerica metodo Runge-Kutta 
 J = [cos(psi(k)) -sin(psi(k)) 0 0;
         sin(psi(k)) cos(psi(k)) 0 0;
         0 0 1 0;
         0 0 0 1];
h_p(:,k) = J*v(:,k+1);

h(:,k+1) = h(:,k)+ UAV_RK4(h(:,k),v(:,k+1),ts);
xu(k+1) = h(1,k+1);
yu(k+1) = h(2,k+1);
zu(k+1) = h(3,k+1);      
psi(k+1) = Angulo(h(4,k));

%% A and B Estimation DMD ONLINE
A_c = (Ae-eye(4))/ts;
B_c = Be/ts;
v_est(:, k+1) = Ae*v_est(:,k)+ Be*vc(:,k);
rho = 1;
if sample >= m
    [Ae,Be,P,G] = DMD_Online(m,v_est,vc,v,P,G,k,rho);
    sample = 0;   
end
sample = sample + 1;



%% 3) Tiempo de máquina   
dt(k) = toc;
A_E(:,k+1) = reshape(A,16,1);
B_E(:,k+1) = reshape(B,16,1);  

A_E_P(:,k+1) = reshape(Ae,16,1);
B_E_P(:,k+1) = reshape(Be,16,1);  
end

save("MIL_test.mat","dt","h","h_p","hd","hd_p","t","v","vc","vcp","vref");
disp('Fin de los calculos')

%*************************************************************************%
%**************ANIMACION SEGUIMIENTO DE TRAYECTORIA **********************%
%% ***********************************************************************%
disp('Animacion RUN')

% 1) Parámetros del cuadro de animacion
figure(1)
axis equal
view(-15,15) % Angulo de vista
cameratoolbar
title ("Simulacion")

% 2) Configura escala y color del UAV
Drone_Parameters(0.02);
H1 = Drone_Plot_3D(xu(1),yu(1),zu(1),0,0,psi(1));hold on


% c) Gráfica de la trayectoria deseada
plot3(xd,yd,zd,'--')
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');

% 5) Simulación de movimiento del manipulador aéreo
for k=1:400:length(t)  
% a) Eliminas los dibujos anteriores del manipulador aéreo
delete(H1);
H1 = Drone_Plot_3D(xu(k),yu(k),zu(k),0,0,psi(k)); hold on
% b) Gráfica la posición deseada vs actual en cada frame
plot3(xu(1:k),yu(1:k),zu(1:k),'r')
hold on
plot3(xd(1:k),yd(1:k),zd(1:k),'b')

pause(0.1)
end

disp('FIN Simulación RUN')  

%%
%******************************************************************************************************************
%********************************************* GR�?FICAS ***********************************************************
%% ****************************************************************************************************************


% 2) Cálculos del Error
figure(2)
hxe= xd - xu;
hye= yd - yu;
hze= zd - zu;
psie= Angulo(psid-psi);
plot(hxe), hold on, grid on
plot(hye)
plot(hze)
plot(psie)
legend("hxe","hye","hze","psie")
title ("Errores de posicion")

% 3) Posiciones deseadas vs posiciones reales del extremo operativo del manipulador aéreo
figure(3)

subplot(4,1,1)
plot(t(1,1:length(xd)),xd)
hold on
plot(t(1,1:length(xu)),xu)
legend("xd","hx")
ylabel('x [m]'); xlabel('s [ms]');
title ("Posiciones deseadas y reales del extremo operativo del manipulador aéreo")

subplot(4,1,2)
plot(t(1,1:length(yd)),yd)
hold on
plot(t(1,1:length(yu)),yu)
legend("yd","hy")
ylabel('y [m]'); xlabel('s [ms]');

subplot(4,1,3)
plot(t(1,1:length(zd)), zd)
hold on
plot(t(1,1:length(zd)), zu)
grid on
legend("zd","hz")
ylabel('z [m]'); xlabel('s [ms]');

subplot(4,1,4)
plot(Angulo(psid))
hold on
plot(psi)
legend("psid","psi")
ylabel('psi [rad]'); xlabel('s [ms]');

% 3) Posiciones deseadas vs posiciones reales del extremo operativo del manipulador aéreo
figure(4)
plot(t(1,1:length(vc)),vc(1,1:end))
hold on
plot(t(1:length(v)),v(1,1:end))
hold on
plot(t(1:length(vref)),vref(1,1:end))
hold on
plot(t(1:length(v_est)),v_est(1,1:end))
legend("ulc","ul","ul_{ref}","ul_{est}")
ylabel('x [m/s]'); xlabel('s [ms]');
title ("Posiciones deseadas y reales del UAV")

figure(5)
plot(t(1:length(vc)),vc(2,1:end))
hold on
plot(t(1:length(v)),v(2,1:end))
hold on
plot(t(1:length(vref)),vref(2,1:end))
hold on
plot(t(1:length(v_est)),v_est(2,1:end))
legend("umc","um","um_{ref}","um_{est}")
ylabel('y [m/s]'); xlabel('s [ms]');

figure(6)
plot(t(1:length(vc)),vc(3,1:end))
hold on
plot(t(1:length(v)),v(3,1:end))
hold on
plot(t(1:length(vref)),vref(3,1:end))
hold on
plot(t(1:length(v_est)),v_est(3,1:end))
legend("unc","un","un_{ref}","un_{est}")
ylabel('z [m/ms]'); xlabel('s [ms]');

figure(7)
plot(t(1:length(vc)),vc(4,1:end))
hold on
plot(t(1:length(v)),v(4,1:end))
hold on
plot(t(1:length(vref)),vref(4,1:end))
hold on
plot(t(1:length(v_est)),v_est(4,1:end))
legend("wc","w","w_{ref}","w_{est}")
ylabel('psi [rad/s]'); xlabel('s [ms]');

figure(8)
subplot(1,2,1)
plot(A_E_P',"-",'linewidth',2)
% legend("wc","w","w_{ref}")
ylabel('psi [rad/s]'); xlabel('s [ms]');
title ("Valores estimados de A")

subplot(1,2,2)
plot(B_E_P',"-",'linewidth',2)
% legend("wc","w","w_{ref}")
ylabel('psi [rad/s]'); xlabel('s [ms]');
title ("Valores estimados de B")

figure(9)
subplot(1,2,1)
plot(A_E',"-",'linewidth',2)
% legend("wc","w","w_{ref}")
ylabel('psi [rad/s]'); xlabel('s [ms]');
title ("Valores perturbados de A")

subplot(1,2,2)
plot(B_E',"-",'linewidth',2)
% legend("wc","w","w_{ref}")
ylabel('psi [rad/s]'); xlabel('s [ms]');
title ("Valores perturbados de B")


% %%%%%%%%%%%%%
%%
figure(10)

figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 4]);
plot(t(1:length(dt)),dt,'Color',[46,188,89]/255,'linewidth',1); hold on
grid on;
legend({'$t_{sample}$'},'Interpreter','latex','FontSize',11,'Orientation','horizontal');
legend('boxoff')
title('$\textrm{Sample Time}$','Interpreter','latex','FontSize',9);
ylabel('$[s]$','Interpreter','latex','FontSize',9);
xlabel('$\textrm{Time}[s]$','Interpreter','latex','FontSize',9);

figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 4]);
plot(t(1:length(Tu)),Tu,'linewidth',1); hold on
grid on;
legend({'$T_{u}$'},'Interpreter','latex','FontSize',11,'Orientation','horizontal');
legend('boxoff')
title('$\textrm{Perturbacion}$','Interpreter','latex','FontSize',9);
ylabel('$[s]$','Interpreter','latex','FontSize',9);
xlabel('$\textrm{Time}[s]$','Interpreter','latex','FontSize',9);

  
figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 4]);
plot(t(1:length(vref)),vref(:,:),'--','linewidth',1); hold on


grid on;
legend({'$T_{u}$'},'Interpreter','latex','FontSize',11,'Orientation','horizontal');
legend('boxoff')
title('$\textrm{Perturbacion}$','Interpreter','latex','FontSize',9);
ylabel('$[s]$','Interpreter','latex','FontSize',9);
xlabel('$\textrm{Time}[s]$','Interpreter','latex','FontSize',9);
