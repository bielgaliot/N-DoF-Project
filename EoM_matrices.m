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

M = [mf 0 0 0 0;
    0 Jf 0 0 0;
    0 0 m1 0 0;
    0 0 0 m2 0;
    0 0 0 0 m3]

C12= c1*(d12+a)+c2*a-c3*(d23-a);
C22= c1*(d12+a)^2+c2*a^2+c3*(d23-a)^2;

C = [c1+c2+c3 C12 -c1 -c2 -c3;
    C12 C22 -c1*(d12+a) -c2*a c3*(d23-a);
    -c1 -c1*(d12+a) c1 0 0;
    -c2 -c2*a 0 c2 0;
    -c3 c3*(d23-a) 0 0 c3]

K12= k1*(d12+a)+k2*a-k3*(d23-a);
K22= k1*(d12+a)^2+k2*a^2+k3*(d23-a)^2;

K = [k1+k2+k3 K12 -k1 -k2 -k3;
    K12 K22 -k1*(d12+a) -k2*a k3*(d23-a);
    -k1 -k1*(d12+a) k1 0 0;
    -k2 -k2*a 0 k2 0;
    -k3 k3*(d23-a) 0 0 k3]

% El codi següent executa el càlcul de l Normal Eigenvalue Problem per determinar les freqüències naturals i els modes de vibració. Aquest fragment no és necessari per a la resta de codis

[psi, lambda] = eig(K,M); %Resolució automàtica del Normal Eigenvalue Problem normalitzada en massa
w_n = sqrt(lambda); %Freqüències naturals en rad/s
f_n = w_n/(2*pi) %Freqüències naturals en Hz
psi %Visualització dels vectors propis, en aquest cas, modes de vibració.