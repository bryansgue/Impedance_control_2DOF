clc, clear, close all;
%% Parametros del Robot y valores fijos
m1= 1.5; m2= 1.5;
l1= 0.4; l2= 0.4; lb= 0.34;
g=9.807;
ts= 0.01; tf= 20;
t= 0:ts:tf;

mul = 10;
% Definicion de Trayectorias
 xd = 0.2 * sin(mul*0.04*t)+0.45;         xd_p = 0.2*mul*0.04*cos(mul*0.04*t);     xd_pp = -0.2*mul*mul*0.04*0.04*sin(mul*0.04*t);
 zd = 0.2 * sin(mul*0.08*t)+0.7;         zd_p = 0.2*mul*0.08*cos(mul*0.08*t);     zd_pp = -0.2*mul*mul*0.08*0.08*sin(mul*0.08*t);               
 

% radio= 0.3;
% xd= radio*sin(t)+0.3;       zd= 0.2*cos(t)+ 0.8;
% xd_p= radio*cos(t);     zd_p= -0.2*sin(t);
% xd_pp= -radio*sin(t);   zd_pp= -0.2*cos(t);

% Condiciones iniciales 
q=[120*pi/180;-0.9];
q_p= [0;0];
k=1;
xEslabon1= l1*cos(q(1,k));
xEslabon2= xEslabon1 + l2*cos(q(1,k)+q(2,k));
zEslabon1= lb + l1*sin(q(1,k));
zEslabon2= zEslabon1 + l2*sin(q(1,k)+q(2,k));
posReal(:,k)= [xEslabon2;zEslabon2];
velReal(:,k)=[0;0];


% Asignacion de constantes para control de fuerza
Ke= diag([7000,0.01]);        % Modulo de elasticidad longitudinal (Mpa,megaPascales:kg/(m*s^2)), Goma: 7, cartilago:24, Tendon:600, Madera:7000 

% Posicion Obstaculo
posObstaculo= [0.51;0];
% posObstaculo= [1;0];

% Definicion de S: Permite seleccionar ejes de control de movimiento y ejes de control de fuerza
S= [0 0;0 1];
Scoma= eye(2)-S;
Kv_x= 20*diag([1,1]); Kp_x= 50*diag([1,1]);
Kv_f= 10*diag([1,1]); Kp_f= 80*diag([1,1]);
fd= 0.25*[100;0.1];
fd_pp= [0;0];
%% Programa general
for k=1:length(t)
    f(:,k)= [0; 0];
%     if posReal(1,k)>posObstaculo(1)            % El obst�culo esta en la posicion 0.5, tiene forma de barra
    if xd(k)>posObstaculo(1)
        f(:,k)= Ke*(posReal(:,k) - posObstaculo);
    else
        f(:,k)= [0; 0];
    end    
    
    % Parametros del modelo din�mico en coordenadas articulares 
    Mq = Matrix_M_SCARA2DOF(m1,m2,l1,l2,q(1,k),q(2,k));
    Cq = Matrix_C_SCARA2DOF(m1,m2,l1,l2,q(1,k),q(2,k),q_p(1,k),q_p(2,k));
    Gq = Matrix_G_SCARA2DOF(m1,m2,l1,l2,q(1,k),q(2,k));
     
   Ja = Jaco_Scara2DOF(l1,l2,q(1,k),q(2,k));
  
    Ja_p= JacoP_Scara2DOF(l1,l2,q(1,k),q(2,k),q_p(1,k),q_p(2,k));
    % Conversion del modelo din�mico a coordenadas cartesianas
    Mcart= inv(Ja)'*Mq*inv(Ja);
    Ccart= inv(Ja)'*Cq*inv(Ja) - inv(Ja)'*Mq*inv(Ja)*Ja_p*inv(Ja);
    Gcart= inv(Ja)'*Gq;

    % Calculo de errores de posicion/velocidad y controlador de movimiento
    err_p(:,k)= [xd_p(k); zd_p(k)] - velReal(:,k);
    err(:,k)= [xd(k); zd(k)] - posReal(:,k);
    uXc(:,k)= S*[xd_pp(:,k); zd_pp(:,k)] + Kv_x*S*err_p(:,k) + Kp_x*S*err(:,k);    % Calculo de u (Ley de control, lazo externo)
    
    % Calculo de errores de fuerza y controlador de fuerza
%     if posReal(1,k)>posObstaculo(1)
    if xd(k)>posObstaculo(1)
        errF(:,k+1)= fd - f(:,k);
        errF_p(:,k+1)= (errF(:,k+1) - errF(:,k))/ts;
        fXc(:,k)= inv(Ke)*(Scoma*fd_pp + Kv_f*Scoma*errF_p(:,k+1) + Kp_f*Scoma*errF(:,k+1));
        S= [0 0;0 1];
    else
        fXc(:,k)=[0;0];
        errF(:,k+1)=[0;0];
        S= eye(2);
    end
%     uc(:,k)= [xd_pp(:,k); zd_pp(:,k)] + Kv_x*err_p(:,k) + Kp_x*err(:,k);    % Calculo de u (Ley de control, lazo externo)
    
    % Suma de controladores
    uc(:,k)= uXc(:,k) + fXc(:,k);
    
    ft(:,k)= Mcart*uc(:,k) + Ccart*velReal(:,k) + Gcart + f(:,k);       % Calculo de ftA (Ley de control, lazo interno)

    % Aplicacion al robot
    accReal(:,k+1)= inv(Mcart)*(ft(:,k) - Ccart*velReal(:,k) - Gcart - f(:,k));
    velReal(:,k+1)= velReal(:,k) + accReal(:,k+1)*ts;
    posReal(:,k+1)= posReal(:,k) + velReal(:,k+1)*ts;
    
    % Obtenci�n de q_p y q para la siguiente iteracion
    q_p(:,k+1)= inv(Ja)*velReal(:,k+1);
    q(:,k+1)= q(:,k) + q_p(:,k+1)*ts;
    
    h1(:,k) = CD1_Scara2DOF(l1,l2,q(1,k),q(2,k));
    h2(:,k) = CD2_Scara2DOF(l1,l2,q(1,k),q(2,k));
    
end
close all

% *************************************************************************
% ************************* ANIMACI�N *************************************
% *************************************************************************
%a)  Configuraci�n de la animaci�n del robot BOSCH-SR800 
     scene = figure;        % new figure
     tam = get(0,'ScreenSize');
  %   set(scene,'position',[10 50 1500 900]); % position and size figure in the screen
     xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]')
      axis([-0.15 1 -0.8 0.8 -0.45 0.6]); % Set axis limits 
      view(-15,15)
     axis equal;            % Set axis aspect ratios
     grid minor;            % Display axes grid lines   
     camlight('headlight'); % Iluminaci�n del robot
     material('dull');
     cameratoolbar          % Control de escena por teclado
     
 %b) Animaci�n de movimiento del ROBOT BOSCH-SR800

     %R1 = robotPlot(q1(1),q2(1),z(1),scaleRobot);hold on
     %D  = discoPlot(hxd(1),hyd(1),hzd(1),'g',scaleRobot);hold on
   %  H  = plot3(hx(1),hy(1),hz(1),'c','LineWidth',4);
     title('ANIMACI�N DE MOVIMIENTO - CONTROLADOR DIN�MICO PD')
P0=[0 0 0];      

%plot3(hxd,hyd,hzd,'*','LineWidth',4);
%line([hxd P0(1)] ,[hyd  P0(2)],[hzd P0(3)],'Color','blue','LineStyle','-','LineWidth',1); hold on
L1 = plot(0,0); hold on
L2 = plot(0,0); hold on
H1  = plot(0,0);
H2  = plot(0,0);
O1  = plot(0,0);
HD  = plot(0,0);
H3  = plot(0,0);
patch([posObstaculo(1) posObstaculo(1) posObstaculo(1)+0.01 posObstaculo(1)+0.01],[0 0.7 0.7 0],'black')
plot(xd(1,60:end),zd(1,60:end)-0.34,'-','LineWidth',0.5,'Color','blue');


Fr = getframe(gcf);
filename = 'Tra_Hibrido_2DOF.gif';
for n = 1:5:length(t)
    drawnow
    axis equal; 
    axis([-0.1 0.8 -0.05 0.7 ]); % Set axis limits 
    view(0,90)
 %   delete(O1)
    delete(L1)
    delete(L2)
     
    delete(H1)   
    delete(H2)
    delete(H3)
   % delete(HD)
%    O1 = patch([posObstaculo(1,n) posObstaculo(1,n) posObstaculo(1,n)+0.01 posObstaculo(1,n)+0.01],[0.35 0.15 0.35 0.15],'black')
    %if n > 10
        H3  = plot(h2(1,2:n),h2(2,2:n),'.','LineWidth',0.1);
    %end
    H1  = plot(h1(1,n),h1(2,n),'*','LineWidth',4,'Color','red');
    H2  = plot(h2(1,n),h2(2,n),'*','LineWidth',4,'Color','blue');
    % H1  = plot3(hx(2:n),hy(2:n),hz(2:n),'--','LineWidth',1);
    L1 = line([h1(1,n) P0(1)] ,[h1(2,n)  P0(2)],'Color','red','LineStyle','-'); hold on
    L2 = line([h2(1,n) h1(1,n)] ,[h2(2,n)  h1(2,n)],'Color','blue','LineStyle','-'); hold on
    
    
        Fr = getframe(1);
    
    im = frame2im(Fr);
    [imind,cm] = rgb2ind(im,256);
    
    if n==1
        imwrite(imind,cm,filename,'gif','Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    
    pause(ts)
end
% *************************************************************************
% *************************** GR�FICAS ************************************
%% ************************************************************************
% figure
%     subplot(2,1,1)
%         plot(err(1,:),'r','lineWidth',2); hold on
%         plot(err(2,:),'b','lineWidth',2); hold on ;
%         grid minor
% %        plot(t,hze,'r','lineWidth',3); 
%         legend('q_1_e','q_2_e','h_z_e')
%         xlabel('Time [s]');ylabel('[ m ]');
%         title('Errores de Control')
%     subplot(2,1,2)
%         plot(xd(1,:) - posReal(1,1:end-1),'m','lineWidth',2); hold on
%         plot(zd(1,:) - posReal(2,1:end-1),'c','lineWidth',2); hold on 
%          
%         grid minor
%         legend('h_x','h_y','h_z','xd','yd','zd')
%         xlabel('Time [s]');
%         ylabel('[m]');
%         title('Evolucion del Extremo Operativo')
        
 %% Graficas
 
 figure (2)

plot(err(1,1:end),'b','lineWidth',2); hold on
plot(err(2,1:end),'r','lineWidth',2); hold on 

legend('error_x','error_y')
xlabel('Sample [k*ts]');
ylabel('[m]');
title('Errores de Control Hibrido')
%%
figure (3)

plot(errF(1,1:end),'b','lineWidth',2); hold on
plot(errF(2,1:end),'r','lineWidth',2); hold on 

legend('error_x','error_y')
xlabel('Sample [k*ts]');
ylabel('[m]');
title('Errores de fuerza Hibrido')
%%
figure (4)
plot(ft(1,:),'b','lineWidth',1); hold on
plot(ft(2,:),'r','lineWidth',1); hold on 

legend('\tau_x','\tau_y')
xlabel('Sample [k*ts]');
ylabel('[m]');
title('Acciones de control')
%%
figure (5)
title('Fuerzas del medio')
plot(f(1,:),'b','lineWidth',1); hold on
plot(f(2,:),'r','lineWidth',1); hold on 

legend('f_x','f_y')
xlabel('Sample [k*ts]');
ylabel('[Nm]');
title('Fuerzas del medio')