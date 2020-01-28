clear all
% Definim constants
mf = 25000;
m1 = 80;
m2 = 80;
m3 = 65;
Jf = 400000;
k1 = 350000;
k2 = 380000;
k3 = 263000;
c1 = 8000;
c2 = 9000;
c3 = 7000;
d12=1.5;
d23=7.5;
a = 1.5;
kc = 580000000;
%Variem Kx i v
kx=0.2
v=30
A1=[]
A2=[]
A3=[]

for kx=0.1:0.005:0.5
%definim equacions
syms zf(t) phi(t) z1(t) z2(t) z3(t)
ode1 =mf*diff(zf,2)+(c1+c2+c3)*diff(zf)+(k1+k2+k3)*zf-c1*diff(z1)-k1*z1-c2*diff(z2)-k2*z2-c3*diff(z3)-k3*z3+(c1*(d12+a)+c2*a-c3*(d23-a))*diff(phi)+(k1*(d12+a)+k2*a-k3*(d23-a))*phi == 0
ode2 =Jf*diff(phi,2)+c1*(d12+a)+c2*a-c3*(d23-a)*diff(zf)+(k1*(d12+a)+k2*a-k3*(d23-a))*zf-c1*(d12+a)*diff(z1)-k1*(d12+a)*z1-c2*a*diff(z2)-k2*a*z2+c3*(d23-a)*diff(z3)+k3*(d23-a)*z3+(c1*(d12+a)^2+c2*a*a+c3*(d23-a)^2)*diff(phi)+(k1*(d12+a)^2+k2*a*a+k3*(d23-a)^2)*phi == 0
ode3 = m1*diff(z1,2)+c1*diff(z1)+k1*z1-k1*zf-c1*diff(zf)-c1*(d12+a)*diff(phi)-k1*(d12+a)*phi ==  -kc*(sin(kx*v*t)+z1)^(3/2)
ode4 = m2*diff(z2,2)+c2*diff(z2)+k2*z2-k2*zf-c2*diff(zf)-c2*(a)*diff(phi)-k2*a*phi == -kc*(sin(kx*(v*t+d12))+z2)^(3/2)
ode5 = m3*diff(z3,2)+c3*diff(z3)+k3*z3-k3*zf-c3*diff(zf)+c1*(d23-a)*diff(phi)+k3*(d23-a)*phi == -kc*(sin(kx*(v*t+d12+d23))+z3)^(3/2)
[V,S] = odeToVectorField([ode1,ode2,ode3,ode4,ode5]);
vars = [zf(t); phi(t); z1(t);z2(t);z3(t)];
M = matlabFunction(V,'vars', {'t','Y'});
%Solucionem l'edo amb les condicions inicials.
sol = ode45(M,[0.01 20],[0.0028 0 -0.0007109 0 0.003 0 0.0024 0 0.2498 0]); %cond inicials
tValues = linspace(0.1,20,5000);
ValorsZf = deval(sol,tValues,9);
ValorsPhi = deval(sol,tValues,3);
ValorsZ1 = deval(sol,tValues,1);
ValorsZ2 = deval(sol,tValues,5);
ValorsZ3 = deval(sol,tValues,7);
Z1mod=max(real(ValorsZ1))
Z2mod=max(real(ValorsZ2))
Z3mod=max(real(ValorsZ3))
ref1=10000;
ref2=10000;
ref3=10000;
for i=1000:3000
    Z1angle=angle(ValorsZ1(i));
    Z2angle=angle(ValorsZ2(i));
    Z3angle=angle(ValorsZ3(i));
    temp1=abs((-0.0028)/(1+Z1mod*exp(1i*Z1angle)));
    temp2=abs((-0.003)/(exp(1i*kx*d12)+Z2mod*exp(1i*Z2angle)));
    temp3=abs((-0.0024)/(exp(1i*kx*(d12+d23))+Z3mod*exp(1i*Z3angle)));
    if temp1<ref1
        ref1=temp1;
    end
    if temp2<ref2
        ref2=temp2;
    end
    if temp3<ref3
        ref3=temp3;
    end
end
A1=[A1 ref1];
A2=[A2 ref2];
A3=[A3 ref3];
end

%Grafiquem la resposta   
figure;
plot([0.1:0.005:0.5],A1,'r')
hold on;
plot([0.1:0.005:0.5],A2,'b')
plot([0.1:0.005:0.5],A3,'g')
xlabel('$k_x$ $(rad/m)$','Interpreter','Latex','Fontsize',12)
ylabel('$A_{r,\,m\grave{a}x}$','Interpreter','Latex','Fontsize',13)
xlim([0.1 0.5]);
title('$A_r$ en  funci\''o de $k_x$','Interpreter','Latex','Fontsize',12)
legend({'$A_r^{(1)}$','$A_r^{(2)}$','$A_r^{(3)}$','$A_{r,\,m \acute{\imath}n}$'},'Interpreter','Latex','Fontsize',12)
grid
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.5;
print(gcf,'lineal1','-dpng','-r1000');
