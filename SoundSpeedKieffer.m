

%calculates sound speed of water/steam mixture in both chemical equilibrium
% and nonequilibrium cases, used to make figure 10a in Karlstrom et al., 2013
%implements Kieffer 1973 for equilibrium (corrects typo in that work) 
%and Kumagai and Chouet for nonequilibrium

%outputs nonequilibrium and equilibirum sound speeds of water + steam
%inputs steam mass fraction and pressure in bars
clear all
close all



%atmospheric pressure in bar
Press=1;

%use Matlab implementation of steam tables to look up
%densities of saturated liquid and vapor as function of pressure
 Rhog = XSteam('rhoV_p',Press);
 Rhol=  XSteam('rhoL_p',Press); 
 

%% case 1: nonequilibrium response for steam

F=logspace(-4,0,500);

Pr = .1e6 ; %reference pressure in Pa
T0=373; %reference temp at 1 bar
R0=462;% %J/K*kg
gamma = 1.25;
%M0=28.98; % for steam
%rho0=0.96; %g/cm^3
%gamma=1.40; %for adiabatic process
K = 2.2*10^6 *1e5; %bulk modulus of water in Pa
Rho_lr=XSteam('rhoL_p',Pr/1e5);

 RhogNE = Press*1e5/(R0*T0); 
 
 RholNE = Rho_lr*exp((Pr-Press*1e5)/K);

 G=T0*R0/(RhogNE^(gamma-1));

Xne = F.*RholNE./(RhogNE+F*RholNE); %get volume fraction from mass fraction


%Rhol=1 ; %ref density of liquid g/cm^3 at 1 bar

for i=1:length(F)

C0(i)= F(i)*Rho_lr*(G/(1e5*Press))^(1/gamma);    
    %nonequilibrium sound speed: there is no exchange of mass between phases on
%timescale of acoustic wave propagation
C_NE(i) = (C0(i)+exp((Pr-Press*1e5)/K))*(sqrt((1+ F(i))*Rho_lr)*sqrt(C0(i)/(gamma*Press*1e5)+...
    1/K *exp((Pr-Press*1e5)/K)))^(-1);
end

%keyboard

%% case 2: equilibrium sound speed


%P=linspace(.5,2,100);


%Ent=logspace(-.001,.86,3000);


Ent=linspace(1,7.5,5000);


for iP=1:length(Ent)
    %Hsat(iP)=XSteam('hL_p',P(iP)); %enthalpy at saturation KJ/Kg
    Hps(iP) = XSteam('h_ps',Press,Ent(iP))*1000; %enthalpy(p,s) J/Kg
    Vlsat(iP)=XSteam('vL_p',Press); %specific liquid volume at saturation m^3/kg
    Vgsat(iP)=XSteam('vV_p',Press); %specific vapor volume at saturation m^3/kg
    Vps(iP) = XSteam('v_ps',Press,Ent(iP)); %specific volume of mixture
    Xv(iP) = XSteam('x_ps',Press,Ent(iP)); %vapor fraction  
    Vv(iP) = XSteam('vx_ps',Press,Ent(iP)); %vapor volume fraction  

    VoidFrac(iP)=Xv(iP)*Rhol/(Rhog+Xv(iP)*Rhol);
 
 %mass fraction vapor
 Mx(iP) = Rhog*Vv(iP)/(Rhog*Vv(iP)+Rhol*(1-Vv(iP)));

L(iP) = LatentHeat(Press)*1000; %latent heat in J/Kg

%form derivatives, use pressure increment of 0.001 bar, centered
%differences
Pinc=.01;
%for iP=2:length(P) 
    dVldP_sat(iP) = (XSteam('vL_p',Press+Pinc)-XSteam('vL_p',Press-Pinc))/(2*Pinc*1e5); 
    %derivative of saturation liquid specific volume
    dVgdP_sat(iP) = (XSteam('vV_p',Press+Pinc)-XSteam('vV_p',Press-Pinc))/(2*Pinc*1e5); 
    %derivative of saturation vapr specific volume 
    dHldP_sat(iP) = 1000*(XSteam('hL_p',Press+Pinc)-XSteam('hL_p',Press-Pinc))/(2*Pinc*1e5); 
    %derivatice of liquid enthalpy at saturation as fctn of pressure
    dLdP(iP) = 1000*(LatentHeat(Press+Pinc)-LatentHeat(Press-Pinc))/(2*Pinc*1e5); 
    %derivative of latent heat as fctn of pressure
    
 
 %molar mass of water
 %Mm = .01802; %kg/mol
 
    
  %find the specific volume of the mixture, assuming linear mixing
 %Vm(iP) = (Xv(iP)/Rhog + (1-Xv(iP))/Rhol);

Vm(iP) = XSteam('v_ps',Press,Ent(iP)); 


%this is the equilibrium sound speed of Kieffer 
C_E(iP) = sqrt(-Vm(iP)^2  / (((1-Xv(iP))*dVldP_sat(iP) + (Xv(iP))* dVgdP_sat(iP) + ...
    (1/Rhog - 1/Rhol)*(Vm(iP)/L(iP) - 1/L(iP) * dHldP_sat(iP) - ...
    (Xv(iP))/L(iP) *dLdP(iP)))));

Corr(iP)=(1/Rhog - 1/Rhol);


%sound speed of Wood 1941 from Lorenz

%C_W(i) = ((fv/vg +(1-fv)/vl)^.5 * (fv*vg/cg^2 + (1-fv)*vl/cl^2)^.5)^-1;


end


figure
loglog(F,C_NE,'b--')
hold on
loglog(Xv(Xv>0&Xv<1),C_E(Xv>0&Xv<1),'b')
xlim([1e-4 1]);
hold off
xlabel('Mass fraction of steam')
ylabel('Sound speed (m/s)')

%plot(Mx,Corr,'.')
%plot(Mx,(Vm./L - 1./L .* dHldP_sat - (Mx)./L .*dLdP))

%close all


