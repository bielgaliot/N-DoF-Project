H = zeros(5,5,1000);
w_vect = logspace(-1,3,1000);
for count=1:1000;
    w=w_vect(count);   
    H(:,:,count)=inv(-w^2*M+w*1*i*C+K);
end

hold off
t = tiledlayout(5,5); %Definim una matriu de gràfics 5x5
for a=1:5
    for b=1:5
        a_char = num2str(a);
        b_char = num2str(b);
        nexttile %Indiquem a que a cada execució del for es representi la següent casella
        loglog(w_vect, squeeze(abs(H(a,b,:))), 'Displayname', strcat(a_char, b_char));
        xlim([10^-1 10^3])
        grid
    end
end

t.Padding = 'none'; %Reduïm l'espaiat entre gràfics
t.TileSpacing = 'none'

hold off