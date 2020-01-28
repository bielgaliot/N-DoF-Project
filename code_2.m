clear all
mf = 2.5E4;
m1 = 80;
m2 = 80;
m3 = 65;
Jf = 4E5;
k1 = 3.5E5;
k2 = 3.8E5;
k3 = 2.63E5;
d12= 1.5;
d23= 7.5;
c1 = 8E3;
c2 = 9E3;
c3 = 7E3;
a = 1.5;
kc = 5.8E8;
g=9.81;
v=30;

M = [mf 0 0 0 0;
    0 Jf 0 0 0;
    0 0 m1 0 0;
    0 0 0 m2 0;
    0 0 0 0 m3];

C12= c1*(d12+a)+c2*a-c3*(d23-a);
C22= c1*(d12+a)^2+c2*a^2+c3*(d23-a)^2;

C = [c1+c2+c3 C12 -c1 -c2 -c3;
    C12 C22 -c1*(d12+a) -c2*a c3*(d23-a);
    -c1 -c1*(d12+a) c1 0 0;
    -c2 -c2*a 0 c2 0;
    -c3 c3*(d23-a) 0 0 c3];

K12= k1*(d12+a)+k2*a-k3*(d23-a);
K22= k1*(d12+a)^2+k2*a^2+k3*(d23-a)^2;

K = [k1+k2+k3 K12 -k1 -k2 -k3;
    K12 K22 -k1*(d12+a) -k2*a k3*(d23-a);
    -k1 -k1*(d12+a) k1 0 0;
    -k2 -k2*a 0 k2 0;
    -k3 k3*(d23-a) 0 0 k3];


%calculem eigenvec i eigenval
Minv=inv(M);
A=Minv*K;
[eigvec,eigval]=eig(K,M);
%per normalitzar els vectors
Matriunorm=[];
for i=1:5
%separem el vector propi columna
Vector=eigvec(:,i);
%convertim la seva primera component a 1 i la resta el que sigui
Vectornorm=(1/(eigvec(1,i))*Vector);
%incorporem a la matriu de mode shapes (harmònics)
Matriunorm=[Matriunorm Vectornorm];
end

%per transformar els eigenvalues a velocitat angular en rad/s
omega_rads=sqrt(diag(eigval));

%per transformar les velocitats angulars a freqüències en Hz
freqs=omega_rads/(2*pi);

%matriu de massa modal
Mtilde=diag(diag((Matriunorm.')*M*Matriunorm));
%matriu elàstica modal
Ktilde=diag(diag((Matriunorm.')*K*Matriunorm));
%matriu d'esmorteïment (damping) modal
Ctilde=(Matriunorm.')*C*Matriunorm;
Ctildeaprox=diag( diag(Ctilde) );
Ctilde2=(eigvec.')*C*eigvec;

novafreq=(1/(2*pi))*diag(sqrt(Ktilde/Mtilde) );


SoE=@(x)[ 
    (k1+k2+k3)*x(1)+(k1*(a+d12)+k2*a-k3*(d23-a))*x(2)-k1*x(3)-k2*x(4) - k3*x(5)-mf*g;
    ((a+d12)*k1+a*k2-(d23-a)*k3)*x(1)+((a+d12)^2*k1+a^2*k2+(d23-a)^2*k3) *x(2)-(a+d12)*k1*x(3)-a*k2*x(4)+(d23-a)*k3*x(5);
    -k1*x(1)-k1*(a+d12)*x(2)+k1*x(3)+kc*x(3)^1.5-m1*g;
    -k2*x(1)-k2*a*x(2)+k2*x(4)+kc*x(4)^1.5-m2*g;
    -k3*x(1)+k3*(d23-a)*x(2)+k3*x(5)+kc*x(5)^1.5-m3*g
];
x = fsolve(SoE,[0; 0; 0; 0; 0]);
x = round(x,7);
zfs=x(1);
tfs=x(2); 
ds_1 = x(3); 
ds_2 = x(4); 
ds_3 = x(5);

%Càlcul linealitzat

M_new=M;
C_new=C;
DeltaMatrix=[0 0 0 0 0;0 0 0 0 0;0 0 ds_1^0.5 0 0;
    0 0 0 ds_2^0.5 0;0 0 0 0 ds_3^0.5];
K_new=K+1.5*kc*DeltaMatrix;

k=0.001:0.001:0.6;
A_max=zeros(length(k),1);
A_r=zeros(length(k),1);
A_r1=zeros(length(k),1);
A_r2=zeros(length(k),1);
A_r3=zeros(length(k),1);
for i=1:length(k)
    w_exc=k(i)*v;
    Z=-w_exc^2*M_new+1i*w_exc*C_new+K_new;
    H=inv(Z);
    F=-1.5*kc*[0;0;ds_1^0.5;ds_2^0.5*exp(1i*k(i)*d12);
        ds_3^0.5*exp(1i*k(i)*(d12+d23))];
    X=H*F;
    A_r1(i)=abs(-ds_1/(1+X(3)));
    A_r2(i)=abs(-ds_2/(exp(1i*k(i)*d12)+X(4)));
    A_r3(i)=abs(-ds_3/(exp(1i*k(i)*(d12+d23))+X(5)));
    A_r(i)=min([A_r1(i) A_r2(i) A_r3(i)]);   
end

% figure;
% plot(k,A_r1,'r')
% hold on;
% plot(k,A_r2,'b')
% plot(k,A_r3,'g')
% plot(k,A_r,'k')
% xlabel('$k_x$ $(rad/m)$','Interpreter','Latex','Fontsize',12)
% ylabel('$A_{r,\,m\grave{a}x}$','Interpreter','Latex','Fontsize',13)
% xlim([0 k(end)]);
% ylim([-1 6])
% title('$A_r$ en  funci\''o de $k_x$','Interpreter','Latex','Fontsize',12)
% legend({'$A_r^{(1)}$','$A_r^{(2)}$','$A_r^{(3)}$','$A_{r,\,m \acute{\imath}n}$'},'Interpreter','Latex','Fontsize',12)
% 
% grid
% grid minor
% 
% ax = gca;
% ax.GridColor = [0, 0, 0];
% ax.GridAlpha=0.3;
% ax.MinorGridColor = [0, 0, 0];
% ax.MinorGridAlpha=0.5;
% print(gcf,'lineal1','-dpng','-r1000');
[A_rmin,index]=min(A_r);

k_calc=k(index);
A_rcal=A_rmin;

w_calc=k_calc*v;
Z_val=-w_calc^2*M_new+1i*w_calc*C_new+K_new;
F_val=-1.5*kc*A_rcal*[0;0;ds_1;ds_2*exp(1i*k_calc*d12);
        ds_3*exp(1i*k_calc*(d12+d23))];
X=Z_val\F_val;

X_norm=abs(X);
X_ang=angle(X);
fin_time=10*pi/w_calc;
t=0:0.001*fin_time:fin_time;
zf=X_norm(1)*sin(w_calc*t*X_ang(1));
tf=X_norm(2)*sin(w_calc*t*X_ang(2));
z1=X_norm(3)*sin(w_calc*t*X_ang(3));
z2=X_norm(4)*sin(w_calc*t*X_ang(4));
z3=X_norm(5)*sin(w_calc*t*X_ang(5));
rmsz=rms(-w_calc^2*zf);
rmst=rms(-w_calc^2*tf);
rmsw1=rms(-w_calc^2*z1);
rmsw2=rms(-w_calc^2*z2);
rmsw3=rms(-w_calc^2*z3);


figure('Renderer', 'painters', 'Position', [10 10 725 1300])
sgtitle(strcat("Resposta temporal del fuselatge per a $k_x=$",num2str(k_calc)," $rad/m$"),'Interpreter','Latex','Fontsize',12);
subplot(5,1,1);
plot(t,-w_calc^2*zf,'b');
hold on
plot(t,ones(length(t),1)*rmsz,'r');
xlim([0 fin_time]);
xlabel('$t$ $[s]$','Interpreter','Latex','Fontsize',12);
ylabel('$\ddot{z}_f$ $[m/s^2]$','Interpreter','Latex','Fontsize',12)
legend('$\ddot{z}_f(t)$',strcat('RMS= ',num2str(rmsz),' $m/s^2$') ,'Interpreter','Latex','Location','Southeast','Fontsize',10)
grid on

ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.5;

subplot(5,1,2);
plot(t,-w_calc^2*tf,'b');
hold on
plot(t,ones(length(t),1)*rmst,'r');
xlim([0 fin_time]);
xlabel('$t$ $[s]$','Interpreter','Latex','Fontsize',12);
ylabel('$\dot{\varphi}_f$ $[s^{-2}]$','Interpreter','Latex','Fontsize',12)
legend('$\ddot{\varphi}_f(t)$',strcat('RMS= ',num2str(rmst),' $s^{-2}$') ,'Interpreter','Latex','Location','Southeast','Fontsize',10)
grid on

ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.5;

subplot(5,1,3);
plot(t,-w_calc^2*z1,'b');
hold on
plot(t,ones(length(t),1)*rmst,'r');
xlim([0 fin_time]);
xlabel('$t$ $[s]$','Interpreter','Latex','Fontsize',12);
ylabel('$\ddot{z}_1$ $[m/s^2]$','Interpreter','Latex','Fontsize',12)
legend('$\ddot{z}_1(t)$',strcat('RMS= ',num2str(rmsw1),' $m/s^2$') ,'Interpreter','Latex','Location','Southeast','Fontsize',10)
grid on

ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.5;

subplot(5,1,4);
plot(t,-w_calc^2*z2,'b');
hold on
plot(t,ones(length(t),1)*rmst,'r');
xlim([0 fin_time]);
xlabel('$t$ $[s]$','Interpreter','Latex','Fontsize',12);
ylabel('$\ddot{z}_2$ $[m/s^2]$','Interpreter','Latex','Fontsize',12)
legend('$\ddot{z}_2(t)$',strcat('RMS= ',num2str(rmsw2),' $m/s^2$') ,'Interpreter','Latex','Location','Southeast','Fontsize',10)
grid on

ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.5;

subplot(5,1,5);
plot(t,-w_calc^2*z3,'b');
hold on
plot(t,ones(length(t),1)*rmst,'r');
xlim([0 fin_time]);
xlabel('$t$ $[s]$','Interpreter','Latex','Fontsize',12);
ylabel('$\ddot{z}_3$ $[m/s^2]$','Interpreter','Latex','Fontsize',12)
legend('$\ddot{z}_3(t)$',strcat('RMS= ',num2str(rmsw3),' $m/s^2$') ,'Interpreter','Latex','Location','Southeast','Fontsize',10)
grid on

ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.5;
%print(gcf,'lineal2','-dpng','-r1200');