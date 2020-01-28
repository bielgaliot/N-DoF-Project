%Es soluciona sistema per trobar la posició d'equilibri
Sistema=@(x)[ 
    (k1+k2+k3)*x(1)+(k1*(a+d12)+k2*a-k3*(d23-a))*x(2)-k1*x(3)-k2*x(4) - k3*x(5)-mf*g;
    ((a+d12)*k1+a*k2-(d23-a)*k3)*x(1)+((a+d12)^2*k1+a^2*k2+(d23-a)^2*k3) *x(2)-(a+d12)*k1*x(3)-a*k2*x(4)+(d23-a)*k3*x(5);
    -k1*x(1)-k1*(a+d12)*x(2)+k1*x(3)+kc*x(3)^1.5-m1*g;
    -k2*x(1)-k2*a*x(2)+k2*x(4)+kc*x(4)^1.5-m2*g;
    -k3*x(1)+k3*(d23-a)*x(2)+k3*x(5)+kc*x(5)^1.5-m3*g
];
x = fsolve(Sistema,[0; 0; 0; 0; 0]);
zfs=x(1);
tfs=x(2); 
ds_1 = x(3); 
ds_2 = x(4); 
ds_3 = x(5);

%Es defineixen les matrius requerides
M_nova=M;
C_nova=C;
MatriuDelta=[0 0 0 0 0;0 0 0 0 0;0 0 ds_1^0.5 0 0;
    0 0 0 ds_2^0.5 0;0 0 0 0 ds_3^0.5];
K_nova=K+1.5*kc*MatriuDelta;

%S'estudia l'efecte de l'evolució de kx
k=0.001:0.001:1;
A_max=zeros(length(k),1);
A_r=zeros(length(k),1);
A_r1=zeros(length(k),1);
A_r2=zeros(length(k),1);
A_r3=zeros(length(k),1);
for i=1:length(k)
    %La freqüència angular d'excitació
    omega=k(i)*v;
    Z=-omega^2*M_nova+1i*omega*C_nova+K_nova;
    H=inv(Z);
    F=-1.5*kc*[0;0;ds_1^0.5;ds_2^0.5*exp(1i*k(i)*d12);
        ds_3^0.5*exp(1i*k(i)*(d12+d23))];
    X=H*F;
    %Calculem els diferents valors de A a partir de l'aproximació lineal
    A_r1(i)=abs(-ds_1/(1+X(3)));
    A_r2(i)=abs(-ds_2/(exp(1i*k(i)*d12)+X(4)));
    A_r3(i)=abs(-ds_3/(exp(1i*k(i)*(d12+d23))+X(5)));
    A_r(i)=min([A_r1(i) A_r2(i) A_r3(i)]);   
end


%Es realitza la representació gràfica
plot(k,A_r1,'r')
hold on;
plot(k,A_r2,'b')
plot(k,A_r3,'g')
plot(k,A_r,'k')
xlabel('$k_x$ $(rad/m)$','Interpreter','Latex','Fontsize',12)
ylabel('$A_{r,\,m\grave{a}x}$','Interpreter','Latex','Fontsize',13)
title('$A_r$ en  funci\''o de $k_x$','Interpreter','Latex','Fontsize',12)
legend({'$A_r^{(1)}$','$A_r^{(2)}$','$A_r^{(3)}$','$A_{r,\,m \acute{\imath}n}$'},'Interpreter','Latex','Fontsize',12)

%Es busca el valor mínim
[A_rmin,index]=min(A_r);
%S'obtenen els paràmetres associats a l'aproximació lineal
k_lin=k(index);
A_lin=A_rmin;
w_lin=k_lin*v;

Z_lin=-w_lin^2*M_nova+1i*w_lin*C_nova+K_nova;
F_lin=-1.5*kc*A_lin*[0;0;ds_1^0.5;ds_2^0.5*exp(1i*k_lin*d12);
        ds_3^0.5*exp(1i*k_lin*(d12+d23))];
X=Z_lin\F_lin;

X_norm=abs(X);
X_ang=angle(X);
fin_time=5;
t=0:0.001*fin_time:fin_time;
%Es calcula l
zf=X_norm(1)*sin(w_lin*t+X_ang(1));
tf=X_norm(2)*sin(w_lin*t+X_ang(2));
z1=X_norm(3)*sin(w_lin*t+X_ang(3));
z2=X_norm(4)*sin(w_lin*t+X_ang(4));
z3=X_norm(5)*sin(w_lin*t+X_ang(5));
rmsz=rms(-w_lin^2*zf);
rmst=rms(-w_lin^2*tf);
rmsw1=rms(-w_lin^2*z1);
rmsw2=rms(-w_lin^2*z2);
rmsw3=rms(-w_lin^2*z3);


sgtitle(strcat("Resposta temporal dels graus de llibertat per a $k_x=$",num2str(k_calc)," $rad/m$"),'Interpreter','Latex','Fontsize',12);
subplot(5,1,1);
plot(t,-w_calc^2*zf,'b');
hold on
plot(t,ones(length(t),1)*rmsz,'r');
xlim([0 fin_time]);
xlabel('$t$ $[s]$','Interpreter','Latex','Fontsize',12);
ylabel('$\ddot{z}_f$ $[m/s^2]$','Interpreter','Latex','Fontsize',12)
legend('$\ddot{z}_f(t)$',strcat('RMS= ',num2str(rmsz),' $m/s^2$') ,'Interpreter','Latex','Location','Southeast','Fontsize',10)

subplot(5,1,2);
plot(t,-w_calc^2*tf,'b');
hold on
plot(t,ones(length(t),1)*rmst,'r');
xlim([0 fin_time]);
xlabel('$t$ $[s]$','Interpreter','Latex','Fontsize',12);
ylabel('$\ddot{\varphi}_f$ $[s^{-2}]$','Interpreter','Latex','Fontsize',12)
legend('$\ddot{\varphi}_f(t)$',strcat('RMS= ',num2str(rmst),' $s^{-2}$') ,'Interpreter','Latex','Location','Southeast','Fontsize',10)
 
subplot(5,1,3);
plot(t,-w_calc^2*z1,'b');
hold on
plot(t,ones(length(t),1)*rmst,'r');
xlim([0 fin_time]);
xlabel('$t$ $[s]$','Interpreter','Latex','Fontsize',12);
ylabel('$\ddot{z}_1$ $[m/s^2]$','Interpreter','Latex','Fontsize',12)
legend('$\ddot{z}_1(t)$',strcat('RMS= ',num2str(rmsw1),' $m/s^2$') ,'Interpreter','Latex','Location','Southeast','Fontsize',10)
 
subplot(5,1,4);
plot(t,-w_calc^2*z2,'b');
hold on
plot(t,ones(length(t),1)*rmst,'r');
xlim([0 fin_time]);
xlabel('$t$ $[s]$','Interpreter','Latex','Fontsize',12);
ylabel('$\ddot{z}_2$ $[m/s^2]$','Interpreter','Latex','Fontsize',12)
legend('$\ddot{z}_2(t)$',strcat('RMS= ',num2str(rmsw2),' $m/s^2$') ,'Interpreter','Latex','Location','Southeast','Fontsize',10)

subplot(5,1,5);
plot(t,-w_calc^2*z3,'b');
hold on
plot(t,ones(length(t),1)*rmst,'r');
xlim([0 fin_time]);
xlabel('$t$ $[s]$','Interpreter','Latex','Fontsize',12);
ylabel('$\ddot{z}_3$ $[m/s^2]$','Interpreter','Latex','Fontsize',12)
legend('$\ddot{z}_3(t)$',strcat('RMS= ',num2str(rmsw3),' $m/s^2$') ,'Interpreter','Latex','Location','Southeast','Fontsize',10)

%Calculem evolució de la força
fc1=1.5*kc*(ds_1)^0.5*(A_lin*(sin(k_lin*v*t)+z1));
fc2=1.5*kc*(ds_2)^0.5*(A_lin*(sin(k_lin*(v*t+d12)+z2)));
fc3=1.5*kc*(ds_3)^0.5*(A_lin*(sin(k_lin*(v*t+d12+d23)+z3)));
rmsfc1=rms(fc1);
rmsfc2=rms(fc2);
rmsfc3=rms(fc3);

sgtitle(strcat("Resposta temporal de les forces de contacte per a $k_x=$",num2str(k_calc)," $rad/m$"),'Interpreter','Latex','Fontsize',12);
subplot(3,1,1);
plot(t,fc1,'b');
hold on
plot(t,ones(length(t),1)*rmsfc1,'r');
xlim([0 fin_time]);
xlabel('$t$ $[s]$','Interpreter','Latex','Fontsize',12);
ylabel('$f^{(1)}$ $[N]$','Interpreter','Latex','Fontsize',12)
legend('$f^{(1)}(t)$',strcat('RMS= ',num2str(rmsfc1),' $N$') ,'Interpreter','Latex','Location','Southeast','Fontsize',10)

subplot(3,1,2);
plot(t,fc2,'b');
hold on
plot(t,ones(length(t),1)*rmsfc2,'r');
xlim([0 fin_time]);
xlabel('$t$ $[s]$','Interpreter','Latex','Fontsize',12);
ylabel('$f^{(2)}$ $[N]$','Interpreter','Latex','Fontsize',12)
legend('$\ddot{z}_f(t)$',strcat('RMS= ',num2str(rmsfc2),' $N$') ,'Interpreter','Latex','Location','Southeast','Fontsize',10)

subplot(3,1,3);
plot(t,fc3,'b');
hold on
plot(t,ones(length(t),1)*rmsfc3,'r');
xlim([0 fin_time]);
xlabel('$t$ $[s]$','Interpreter','Latex','Fontsize',12);
ylabel('$f^{(3)}$ $[N]$','Interpreter','Latex','Fontsize',12)
legend('$\ddot{z}_f(t)$',strcat('RMS= ',num2str(rmsfc2),' $N$') ,'Interpreter','Latex','Location','Southeast','Fontsize',10)

  