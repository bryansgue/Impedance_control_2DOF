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

% Asignacion de constantes para control de fuerza
Ke= diag([60;0]);        % Modulo de elasticidad longitudinal (Mpa,megaPascales:kg/(m*s^2)), Goma: 7, cartilago:24, Madera:7000 
I= 0.25*diag([1,1]);     % Provee comportamiento suave del extremo ante fuerzas de contacto
D= 5.0*diag([1,1]);     % Valores alto: alta disipacion de energia
K= 18*diag([10,1]);    % Valores altos en direcciones que requieren alta precision de posicion. Bajos en direcciones que req. pequeñas fuerzas de interaccion

k=1;
xEslabon1= l1*cos(q(1,k));
xEslabon2= xEslabon1 + l2*cos(q(1,k)+q(2,k));
zEslabon1= lb + l1*sin(q(1,k));
zEslabon2= zEslabon1 + l2*sin(q(1,k)+q(2,k));

posReal(:,k)= [xEslabon2;zEslabon2];
velReal(:,k)=[0;0];
%% Programa general
for k=1:length(t)

%     if posReal(1,k)>0.5             % El obstáculo esta en la posicion 0.5, tiene forma de barra
%         f(:,k)= Ke*([0.5;0] - posReal(:,k));
%     else
%         f(:,k)= [0; 0];
%     end
    f(:,k)= [0; 0];
    
    % Parametros del modelo dinámico en coordenadas articulares 
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
    
    % Conversion del modelo dinámico a coordenadas cartesianas
    Mcart= inv(Ja)'*Mq*inv(Ja);
    Ccart= inv(Ja)'*Cq*inv(Ja) - inv(Ja)'*Mq*inv(Ja)*Ja*inv(Ja);
    Gcart= inv(Ja)'*Gq;

    % Calculo de errores
    err_p(:,k)= [xd_p(k); zd_p(k)] - velReal(:,k);
    err(:,k)= [xd(k); zd(k)] - posReal(:,k);

    % Calculo del controlador de fuerza
    u(:,k)= inv(Ja)*inv(I)*(I*[xd_pp(:,k); zd_pp(:,k)] + D*err_p(:,k) + K*err(:,k) - I*Ja_p*q_p(:,k) - f(:,k));    % Calculo de u (Ley de control, lazo externo)
    torque(:,k)= Mq*u(:,k) + Cq*q_p(:,k) + Gq + Ja'*f(:,k);       % Calculo de ftA (Ley de control, lazo interno)

    % Aplicacion al robot
    q_pp(:,k+1)= inv(Mq)*(torque(:,k) - Cq*q_p(:,k) - Gq - Ja'*f(:,k));
    
    % Obtención de q_p y q para la siguiente iteracion
    q_p(:,k+1)= q_p(:,k) + q_pp(:,k+1)*ts;
    q(:,k+1)= q(:,k) + q_p(:,k+1)*ts;
    
    velReal(:,k+1)= Ja*q_p(:,k+1);
%     posReal(:,k+1)= posReal(:,k) + velReal(:,k+1)*ts;
    
    xEslabon1= l1*cos(q(1,k));
    xEslabon2= xEslabon1 + l2*cos(q(1,k)+q(2,k));
    zEslabon1= lb + l1*sin(q(1,k));
    zEslabon2= zEslabon1 + l2*sin(q(1,k)+q(2,k));
    posReal(:,k+1)=[xEslabon2; zEslabon2];
end

plot(err(1:300)')


%% Animacion
% close all
qr= q;
% Calculo de las posiciones de los eslabores
xEslabon1= l1*cos(qr(1,:));
xEslabon2= xEslabon1 + l2*cos(qr(1,:)+qr(2,:));
zEslabon1= lb + l1*sin(qr(1,:));
zEslabon2= zEslabon1 + l2*sin(qr(1,:)+qr(2,:));
slab1=[]; slab2=[]; slab=[]; slab3=[]; slab4=[];

figure(11)
hold on, grid on
axis([-1.2 1.2 0 1.5])
plot(xd,zd,':','linewidth',1)
% Ploteo de la base
plot([0 0],[0 lb],'color',[253 106 0]/255,'linewidth',6)
plot(0,lb,'ok','linewidth',6)
% Animacion de lo que se ejecuta
for k= 1:length(xEslabon1)-1
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
% figure(9),
% subplot(3,1,1)
% hold on,grid on
% plot(t,err_q','linewidth',1.5)
% legend('err_{q1}','err_{q2}','Orientation','horizontal');
% set(gca,'FontWeight','Bold')
% title('Error de posicion del E.O.')
% set(gca,'FontWeight','Bold')
% title('Error de posiciones articulares')
% % axis([0 10 -1 1])
% subplot(3,1,2)
% hold on,grid on
% plot(t,err_q_p','linewidth',1.5)
% legend('err_{q1_p}','err_{q2_p}','Orientation','horizontal');
% set(gca,'FontWeight','Bold')
% title('Error de velocidades articulares')
% subplot(3,1,3)
% hold on,grid on
% plot(t,err_q_pp','linewidth',1.5)
% legend('err_{q1_{pp}}','err_{q2_{pp}}','Orientation','horizontal');
% set(gca,'FontWeight','Bold')
% title('Error de aceleraciones articulares')
% xlabel('Tiempo [s]')