clc, clear, close all;
%% Parametros del Robot y valores fijos
m1= 1; m2= 1;
l1= 0.4; l2= 0.526; lb= 0.34;
g=9.807;
ts= 1/100; 
t= 0:ts:10;

% Definicion de Trayectorias
radio= 0.3;
xd= radio*sin(t)+0.3;   zd= 0.2*cos(t)+ 0.8;
xd_p= radio*cos(t);     zd_p= -0.2*sin(t);
xd_pp= -radio*sin(t);   zd_pp= -0.2*cos(t);

xd = 0.5*ones(1,length(t));
zd = 0.25*ones(1,length(t));
xd_p = 0*ones(length(t));
zd_p= 0*ones(length(t));
xd_pp= 0*ones(length(t));
zd_pp= 0*ones(length(t));

% Condiciones iniciales 
q=[0*pi/180;150*pi/180];
%q=[-60*pi/180;-90*pi/180];
q_p= [0;0];

% Asignacion de constantes para control de fuerza
Ke= diag([3000;0]);        % Modulo de elasticidad longitudinal (Mpa,megaPascales:kg/(m*s^2)), Goma: 7, cartilago:24, Madera:7000 
I= 1*diag([1,1]);     % Provee comportamiento suave del extremo ante fuerzas de contacto
D= 7*diag([1,1]);     % Valores alto: alta disipacion de energia
K= 5*diag([10,1]);    % Valores altos en direcciones que requieren alta precision de posicion. Bajos en direcciones que req. peque�as fuerzas de interaccion

I= 3*diag([1,1]);     % Provee comportamiento suave del extremo ante fuerzas de contacto
D= 7*diag([1,1]);     % Valores alto: alta disipacion de energia
K= 5*diag([10,1]);

k=1;

h1(:,k) = CD1_Scara2DOF(l1,l2,q(1,k),q(2,k))  
h2(:,k) = CD2_Scara2DOF(l1,l2,q(1,k),q(2,k));

posReal(:,k)= [h2(1,k);h2(2,k)];
velReal(:,k)=[0;0];

obx = 0.2*sin(2*t).*cos(t)+0.4;
oby = 0*ones(length(t));
% Posicion Obstaculo
posObstaculo(:,1)= [obx(1);oby(1)];

%% Programa general
for k=1:length(t)

    % Calculo de errores
    err_p(:,k)= [xd_p(k); zd_p(k)] - velReal(:,k);
    err(:,k)= [xd(k); zd(k)] - posReal(:,k);
    
    posObstaculo(:,k) = [obx(k);oby(k)];
    
    f(:,k)= [0; 0];
    if posReal(1,k)>posObstaculo(:,k)            % El obst�culo esta en la posicion 0.5, tiene forma de barra
        f(:,k)= D*velReal(:,k)+Ke*(posReal(:,k) - posObstaculo(:,k));
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

    % Calculo del controlador de fuerza
    u(:,k)= [xd_pp(k); zd_pp(k)] + inv(I)*(D*err_p(:,k) + K*err(:,k) - f(:,k));    % Calculo de u (Ley de control, lazo externo)
    ft(:,k)= Mcart*u(:,k) + Ccart*velReal(:,k) + Gcart + f(:,k);       % Calculo de ftA (Ley de control, lazo interno)
    
     %ft(:,k) = [0;0];
    
    ff1 = [1.0256*sign(velReal(1,k));...
        1.7842*sign(velReal(2,k))];
    
%     % Aplicacion al robot;
    accReal(:,k+1)= inv(Mcart)*(ft(:,k) - Ccart*velReal(:,k) - Gcart - 0.0*ff1 - f(:,k));
    velReal(:,k+1) = velReal(:,k) + x_pp_RK4_2DOF(velReal(:,k),Mcart, Ccart, Gcart,ft(:,k), ff1, f(:,k),ts);
    posReal(:,k+1) = posReal(:,k) + x_p_RK4_2DOF(velReal(:,k+1),ts);

    % Obtenci�n de q_p y q para la siguiente iteracion
    q_p(:,k+1)= inv(Ja)*velReal(:,k+1);
     q(:,k+1)= q(:,k) + q_p_RK4_2DOF(velReal(:,k+1),Ja,ts);
   % q(:,k+1) = ts*q_p(:,k)+q(:,k);

%     T = [ cos(q(1,k))*cos(q(2,k)) - sin(q(1,k))*sin(q(2,k)), - cos(q(1,k))*sin(q(2,k)) - cos(q(1,k))*sin(q(1,k)) ;
%          cos(q(1,k))*sin(q(2,k)) + cos(q(2,k))*sin(q(1,k)),   cos(q(1,k))*cos(q(2,k)) - sin(q(1,k))*sin(q(1,k))];
%     ff = [sign(q_p(1,k));...
%           sign(q_p(2,k))];
%      
%     T_ref(:,k) = Ja'*ft(:,k);
%     tf(:,k) = Ja'*f(:,k);
% 
%     qpp = inv(Mq)*(T_ref(:,k)-Cq*q_p(:,k)-Gq-0.1*ff-tf(:,k)); %Aceleraci�n articular de salida
%    
%     q_p(:,k+1) = ts*qpp+q_p(:,k);
%    
%     q(:,k+1) = ts*q_p(:,k)+q(:,k);
%     
    h1(:,k) = CD1_Scara2DOF(l1,l2,q(1,k),q(2,k));
    h2(:,k) = CD2_Scara2DOF(l1,l2,q(1,k),q(2,k));
%     posReal(:,k+1) = h2(:,k); 
       
end
 
%% 
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
     scaleRobot = 1;
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

Fr = getframe(gcf);
filename = 'Pos_Impedance_2DOF.gif';
for n = 1:5:length(t)
    drawnow
    axis equal; 
    axis([-0.1 0.8 -0.4 0.5 -1 1]); % Set axis limits 
    view(0,90)
    delete(O1)
    delete(L1)
    delete(L2)
     
    delete(H1)   
    delete(H2)
    delete(H3)
    delete(HD)
    O1 = patch([posObstaculo(1,n) posObstaculo(1,n) posObstaculo(1,n)+0.01 posObstaculo(1,n)+0.01],[0.35 0.15 0.15 0.35],'black')

    H1  = plot(h1(1,n),h1(2,n),'*','LineWidth',4,'Color','red');
    H2  = plot(h2(1,n),h2(2,n),'*','LineWidth',4,'Color','blue');
    if n > 10
        H3  = plot(h2(1,n-10:n),h2(2,n-10:n),'.','LineWidth',0.1);
    end
    
    L1 = line([h1(1,n) P0(1)] ,[h1(2,n)  P0(2)],'Color','red','LineStyle','-'); hold on
    L2 = line([h2(1,n) h1(1,n)] ,[h2(2,n)  h1(2,n)],'Color','blue','LineStyle','-'); hold on
    HD= plot(xd(k),zd(k),'o','linewidth',1.5);
    xlabel('X [m]');ylabel('Y [ m ]');
    
    Fr = getframe(1);
    
    im = frame2im(Fr);
    [imind,cm] = rgb2ind(im,256);
    
%     if n==1
%         imwrite(imind,cm,filename,'gif','Loopcount',inf);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append');
%     end

    pause(ts)
    
    
    
end
% *************************************************************************
% *************************** GR�FICAS ************************************
%% ************************************************************************
figure (2)

plot(xd(1,:) - posReal(1,1:end-1),'b','lineWidth',2); hold on
plot(zd(1,:) - posReal(2,1:end-1),'r','lineWidth',2); hold on 

legend('error_x','error_y')
xlabel('Time [s]');
ylabel('[m]');
title('Errores de Control')

figure (3)
plot(u(1,:),'b','lineWidth',1); hold on
plot(u(2,:),'r','lineWidth',1); hold on 

legend('u_x','u_y')
xlabel('Time [s]');
ylabel('[m]');
title('Lazo externo')

figure (4)
title('Fuerzas del medio')
plot(f(1,:),'b','lineWidth',1); hold on
plot(f(2,:),'r','lineWidth',1); hold on 

legend('f_x','f_y')
xlabel('Time [s]');
ylabel('[Nm]');
title('Fuerzas del medio')

figure (5)
plot((ft(1,15:end)),'b','lineWidth',1); hold on
plot((ft(2,15:end)),'r','lineWidth',1); hold on 
legend('\tau_x','\tau_y')
xlabel('Time [s]');
ylabel('[Nm]');

title('Fuerzas de control')
 
 
        
      
    
    