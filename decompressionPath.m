%Calculate internal energy
clear all
close all

pmin = 1.0;
pmax = 16;
Pres = 900; %bar
Tres = 205; %celsius
Ttop=100;

P = linspace(pmin,pmax,100);
 P = linspace((pmin),(pmax),500);
T = zeros(length(P),1);
T(length(T)) = Tres;

resEnthalpy = XSteam('hL_T',Tres);%XSteam('h_pT',Pres,Tres);
resEntropy = XSteam('sL_T',Tres);%XSteam('s_pT',Pres,Tres);

topEnthalpy = XSteam('h_pT',pmin,Ttop);
topEntropy = XSteam('s_pT',pmin,Ttop);


T = zeros(length(P),1);
for iP = 1:length(P)
    T(iP) = XSteam('T_ph',P(iP),resEnthalpy);
    SSs(iP) = XSteam('w_ps',P(iP),resEntropy);
    SSH(iP) = XSteam('w_ph',P(iP),resEnthalpy);
    SS2(iP) = XSteam('wV_p',P(iP)); %saturated vapor sound speed
    SS3(iP) = XSteam('wL_p',P(iP)); %saturated liquid sound speed
    fH(iP) = XSteam('x_ph',P(iP),resEnthalpy);
    fs(iP) = XSteam('x_ps',P(iP),resEntropy);
    Vg(iP) = XSteam('vV_p',P(iP));
    Vl(iP) = XSteam('vL_p',P(iP));
    Vph(iP) = XSteam('v_ph',P(iP),resEnthalpy);
    %calculate sound speed
    
end

figure(1)
clf;

hold on
plot((P),T);
xlabel('Pressure (bar)')
%ylabel('T in C')
plot((P),SSH,'r')
plot((P),SSs,'m')
plot((P),SS2,'g')
plot((P),SS3,'k')
figure
%plot((P),Vg)
hold on
plot(P,1./Vl,'r');
%plot(P,Vph,'g')
ylabel('mass fraction steam')
