H = zeros(5,5,1000); %Pre-allocate la matriu H per a totes les freqüències analitzades
w_vect = logspace(-1,3,1000); %Rang de freqüències en distribució logarítmica per optimitzar la representació logarítmica

for count=1:1000;
    w=w_vect(count);   
    H(:,:,count)=inv(-w^2*M+w*1*i*C+K);  %Càlcul de la matriu per a cada freqüència
end

hold off
for a=1:5
    for b=1:5 %Recorrem al matriu casella per casella
        a_char = num2str(a);
        b_char = num2str(b); %índex de cada casella per a la llegenda
        loglog(w_vect, squeeze(abs(H(a,b,:))), 'Displayname', strcat(a_char, b_char));
        hold on
    end
end
grid

%Donem format al gràfic
xlabel('Freq\"u\`encia d''excitaci\''o $\omega$ $(rad/s)$','interpreter','latex')
ylabel('M\`odul de la recept\`ancia $ H_{ij} $', 'interpreter', 'latex')
legend('Location', 'Best')

hold off