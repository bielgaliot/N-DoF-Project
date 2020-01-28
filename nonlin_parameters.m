
kx=0.22;
Ar=0.00182;
%Es defineixen les equacions
syms zf(t) phi(t) z1(t) z2(t) z3(t)
ode1 =mf*diff(zf,2)+(c1+c2+c3)*diff(zf)+(k1+k2+k3)*zf-c1*diff(z1)-k1*z1-c2*diff(z2)-k2*z2-c3*diff(z3)-k3*z3+(c1*(d12+a)+c2*a-c3*(d23-a))*diff(phi)+(k1*(d12+a)+k2*a-k3*(d23-a))*phi == 0
ode2 =Jf*diff(phi,2)+c1*(d12+a)+c2*a-c3*(d23-a)*diff(zf)+(k1*(d12+a)+k2*a-k3*(d23-a))*zf-c1*(d12+a)*diff(z1)-k1*(d12+a)*z1-c2*a*diff(z2)-k2*a*z2+c3*(d23-a)*diff(z3)+k3*(d23-a)*z3+(c1*(d12+a)^2+c2*a*a+c3*(d23-a)^2)*diff(phi)+(k1*(d12+a)^2+k2*a*a+k3*(d23-a)^2)*phi == 0
ode3 = m1*diff(z1,2)+c1*diff(z1)+k1*z1-k1*zf-c1*diff(zf)-c1*(d12+a)*diff(phi)-k1*(d12+a)*phi ==  -Ar*kc*(sin(kx*v*t)+z1)^(3/2)
ode4 = m2*diff(z2,2)+c2*diff(z2)+k2*z2-k2*zf-c2*diff(zf)-c2*(a)*diff(phi)-k2*a*phi == -Ar*kc*(sin(kx*(v*t+d12))+z2)^(3/2)
ode5 = m3*diff(z3,2)+c3*diff(z3)+k3*z3-k3*zf-c3*diff(zf)+c1*(d23-a)*diff(phi)+k3*(d23-a)*phi == -Ar*kc*(sin(kx*(v*t+d12+d23))+z3)^(3/2)
%Es resolen les equacions diferencials
[V,S] = odeToVectorField([ode1,ode2,ode3,ode4,ode5]);
vars = [zf(t); phi(t); z1(t);z2(t);z3(t)];
M = matlabFunction(V,'vars', {'t','Y'});
sol = ode45(M,[0.01 20],[0.0028 0 -0.0007109 0 0.003 0 0.0024 0 0.2498 0]); %cond inicials
tValues = linspace(0.1,20,5000);
ValorsZf = deval(sol,tValues,9);
ValorsPhi = deval(sol,tValues,3);
ValorsZ1 = deval(sol,tValues,1);
ValorsZ2 = deval(sol,tValues,5);
ValorsZ3 = deval(sol,tValues,7);

%Es retalla el tram estacionari del senyal
zf=real(ValorsZf(end/2:end));
tf=real(ValorsPhi(end/2:end));
z1=real(ValorsZ1(end/2:end));
z2=real(ValorsZ2(end/2:end));
z3=real(ValorsZ3(end/2:end));
time=0:0.004:10; %Vector temps

%Es calculen paràmetres d'interés
omega=0.22*30;
rmsz=rms(-omega^2*zf);
rmst=rms(-omega^2*tf);
rmsw1=rms(-omega^2*z1);
rmsw2=rms(-omega^2*z2);
rmsw3=rms(-omega^2*z3);

%El mateix per a la força
fc1=-Ar*kc*(sin(kx*v*time)+z1).^(3/2);
fc2=-Ar*kc*(sin(kx*(v*time+d12))+z2).^(3/2);
fc3=-Ar*kc*(sin(kx*(v*time+d12+d23))+z3).^(3/2);
rmsfc1=rms(fc1);
rmsfc2=rms(fc2);
rmsfc3=rms(fc3);


