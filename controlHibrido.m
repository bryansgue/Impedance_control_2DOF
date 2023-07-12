clc, clear, close all;
%% Parametros del Robot y valores fijos
m1= 5.067; m2= 3.152;
l1= 0.4; l2= 0.526; lb= 0.34;
g=9.807;
ts= 0.01; tf= 10;
t= 0:ts:tf;

% Definicion de Trayectorias
radio= 0.3;
xd= radio*sin(t)+0.3;       zd= 0.2*cos(t)+ 0.8;
xd_p= radio*cos(t);     zd_p= -0.2*sin(t);
xd_pp= -radio*sin(t);   zd_pp= -0.2*cos(t);

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
Kv_x= 20*diag([1,1]); Kp_x= 100*diag([1,1]);
Kv_f= 10*diag([1,1]); Kp_f= 80*diag([1,1]);
fd= [300;0.1];
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
    Mq= [(m1+m2)*l1^2 + m2*l2^2 + 2*m2*l1*l2*cos(q(2,k)), m2*l2^2+m2*l1*l2*cos(q(2,k));
        m2*l2^2+m2*l1*l2*cos(q(2,k))                   , m2*l2^2];
    Cq= [-m2*l1*l2*sin(q(2,k))*q_p(2,k),  -m2*l1*l2*sin(q(2,k))*(q_p(1,k)+q_p(2,k));
        m2*l1*l2*sin(q(2,k))*q_p(1,k) , 0];
    Gq= [(m1+m2)*g*l1*cos(q(1,k)) + m2*g*l2*cos(q(1,k)+q(2,k));
        m2*g*l2*cos(q(1,k)+q(2,k)) ];
    
    % Jacobianita
    Ja=[-l1*sin(q(1,k))-l2*sin(q(1,k)+q(2,k))  -l2*sin(q(1,k)+q(2,k));
    l1*cos(q(1,k))+l2*cos(q(1,k)+q(2,k))    l2*cos(q(1,k)+q(2,k))];

    % Derivada de la Jacobianita (temporal)
    Ja_p= [-l1*cos(q(1,k))*q_p(1,k)-l2*cos(q(1,k)+q(2,k))*(q_p(1,k)+q_p(2,k)), -l2*cos(q(1,k)+q(2,k))*(q_p(1,k)+q_p(2,k));...
    -l1*sin(q(1,k))*q_p(1,k)-l2*sin(q(1,k)+q(2,k))*(q_p(1,k)+q_p(2,k)),        -l2*sin(q(1,k)+q(2,k))*(q_p(1,k)+q_p(2,k))];
    
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
end

%% Animacion
close all
qr= q;
% Calculo de las posiciones de los eslabores
xEslabon1= l1*cos(qr(1,:));
xEslabon2= xEslabon1 + l2*cos(qr(1,:)+qr(2,:));
zEslabon1= lb + l1*sin(qr(1,:));
zEslabon2= zEslabon1 + l2*sin(qr(1,:)+qr(2,:));
slab1=[]; slab2=[]; slab=[]; slab3=[]; slab4=[];

figure(11)
hold on, grid on
patch([posObstaculo(1) posObstaculo(1) posObstaculo(1)+0.1 posObstaculo(1)+0.1],[0 1.5 1.5 0],'black')
axis([-1.2 1.2 0 1.5])
plot(xd,zd,'r:','linewidth',1)
% Ploteo de la base
plot([0 0],[0 lb],'color',[253 106 0]/255,'linewidth',6)
plot(0,lb,'ok','linewidth',6)
% Animacion de lo que se ejecuta
for k= 1:1:length(xEslabon1)-1
    delete(slab1);delete(slab2);delete(slab);delete(slab3); delete(slab4)
    slab1= plot([0 xEslabon1(k)],[lb zEslabon1(k)],'color',[253 106 0]/255,'linewidth',5);
    slab2= plot([xEslabon1(k) xEslabon2(k)],[zEslabon1(k) zEslabon2(k)],'color',[253 106 0]/255,'linewidth',4);
    slab= plot([xEslabon1(k) xEslabon2(k)],[zEslabon1(k) zEslabon2(k)],'ok','linewidth',4);
    slab3= plot([xEslabon2(1:k)],[zEslabon2(1:k)],'m','linewidth',1.5);
    slab4= plot(xd(k),zd(k),'og','linewidth',1.5);
    pause(0.02)
end
set(gca,'FontWeight','Bold')
xlabel('X [m]')
ylabel('Z [m]')


%% Graficacion de Resultados
figure(9),
subplot(2,1,1)
hold on,grid on
plot(t,err','linewidth',1.5)
legend('err_x','err_y');
set(gca,'FontWeight','Bold')
title('Error de posicion del E.O.')

subplot(2,1,2)
hold on,grid on
plot(t,errF(:,1:end-1)','linewidth',1.5)
legend('err_{f_x}','err_{f_y}');
set(gca,'FontWeight','Bold')
title('Error de Fuerzas')

uc(:,k)= uXc(:,k) + fXc(:,k);
figure(10)
subplot(3,1,1)
hold on, grid on
plot(t,uXc','linewidth',1.5)
legend('uP_{x}','uP_{z}');
title('Controlador de Posici�n')
subplot(3,1,2)
hold on, grid on
plot(t,fXc','linewidth',1.5)
legend('uF_{x}','uF_{z}');
title('Controlador de Fuerza')
subplot(3,1,3)
hold on, grid on
plot(t,uc','linewidth',1.5)
legend('u_{x}','u_{z}');
title('Suma de controladores')