clear all
close all

%--------------------------------------------------------
%          Modes of a String Terminated with
%                    EB Beam
%       (String is simply-supported at x=0)
%                Dr M Ducceschi
%             University of Bologna
%                  14 Jul 2022
%--------------------------------------------------------

%-- string parameters
T0    = 147.7 ;            %-- tension [N]
r     = 4.55e-4 ;          %-- radius [m]
E     = 2.5e10 ;           %-- Young's mod [Pa]
rho   = 5535 ;             %-- density [kg/m^3]
L     = 0.69 ;             %-- length [m]

%-- bridge parameters
Eb    = 2e11 ;             %-- Young's mod [Pa]
rhob  = 8e3 ;              %-- density [N/m]
rb    = 1e-3 ;           %-- radius [m]
Lb    = 7e-2 ;              %-- length [m]
zc    = 0.13 ;             %-- connection location [frac: [0-1] ] ;

MaxFreq = 14e3 ;           %-- largest frequency for NR

%--------------------------------------------------------

%-- derived paramters
A     = pi*r^2 ;
I     = pi*r^4/4 ;
EI    = E*I ;
rA    = rho*A ;

Ab     = pi*rb^2 ;
Ib     = pi*rb^4/4 ;
EIb    = Eb*Ib ;
rAb    = rhob*Ab ;

%-- compute beam modes up to MaxOm
MaxOm     = 2*pi*MaxFreq ;
Om        = 0 ;
m         = 0 ;
Wzc       = [] ;
S         = [] ;
while Om <= MaxOm
    m = m + 1 ;
    Om = sqrt(EIb/rAb) * (m^2 * pi^2 / Lb^2) ;
    Wzc = [Wzc; sqrt(2/Lb)*sin(m*pi*zc)] ;
    S   = [S; m^4*pi^4/Lb^4] ;
end

Mmax   = m ;

S = diag(S) ;

% Om = Om(1) ;
% Wzc = Wzc(1) ;
% S = S(1) ;



%--------------------------------------------------------

%-- newton-raphson routine for eigenvalues

om    = 20 * 2 * pi : MaxOm ;
g     = gdef(om,rAb,EIb,Wzc,S) ;
gp    = gpdef(om,rAb,EIb,Wzc,S) ;
Ps    = Psdef(om,T0, EI, rA, L) ;
Psp   = Pspdef(om,T0, EI, rA, L) ;
Psppp = Pspppdef(om,T0, EI, rA, L) ;

%-- eigenvale equation
fom   =  Ps - g.* (-T0*Psp + EI*Psppp) ;
plot(om,fom) ; hold on ; grid on ;
mm = max(abs(fom)) ;
ylim([-10,10])



%-- search for zeros (approximate search)
zom = [] ;
for n = 1 : length(om) - 1

    tt = fom(n) * fom(n+1) ;

    if tt <= 0

        zom = [zom;om(n)] ;

    end

end


%-- initialise zero vector
zerF = zeros(length(zom),1) ;

%-- newton Raphson
for n = 1 : length(zom)

    om = zom(n) ;
    tol = 1 ;
    iter = 0 ;

    while tol > 1e-10 && iter < 20

        iter  = iter + 1 ;

        g        = gdef(om,rAb,EIb,Wzc,S) ;
        gp       = gpdef(om,rAb,EIb,Wzc,S) ;
        Ps       = Psdef(om,T0, EI, rA, L) ;
        Psp      = Pspdef(om,T0, EI, rA, L) ;
        Psppp    = Pspppdef(om,T0, EI, rA, L) ;
        dPs      = dPsdef(om,T0, EI, rA, L) ;
        dPsp     = dPspdef(om,T0, EI, rA, L) ;
        dPsppp   = dPspppdef(om,T0, EI, rA, L) ;


        fom     = Ps - g.*(-T0*Psp + EI*Psppp) ;
        dfom    = dPs - gp.*(-T0*Psp + EI*Psppp) - g.*(-T0*dPsp + EI*dPsppp) ;
        omk     = om - fom/dfom ;

        tol     = abs(omk-om) ;
        om      = omk ;

    end

    g        = gdef(om,rAb,EIb,Wzc,S) ;
    Ps       = Psdef(om,T0, EI, rA, L) ;
    Psp      = Pspdef(om,T0, EI, rA, L) ;
    Psppp    = Pspppdef(om,T0, EI, rA, L) ;
    fom     = Ps - g.*(-T0*Psp + EI*Psppp) ;
    zerF(n) = [0] ;
    if abs(fom) < 1
        zerF(n) = om ;


        line([om,om],[-mm,mm],'linestyle','--','color','k') ;

    end

end



function y = gdef(om,rAb,EIb,Wzc,S)
y     = zeros(1,length(om)) ;
sS    = size(S) ;
Id    = eye(sS(1)) ;
for m = 1 : length(om)
    y(m)     = Wzc.' * ( ( - rAb * om(m).^2 * Id + EIb * S)^(-1) ) * Wzc ;
end
end

function y = gpdef(om,rAb,EIb,Wzc,S)
y     = zeros(1,length(om)) ;
sS    = size(S) ;
Id    = eye(sS(1)) ;
for m = 1 : length(om)
    temp    = - ( - rAb * om(m).^2 * Id + EIb * S)^(-1) *  ( - 2 * rAb * om(m) * Id ) * ( - rAb * om(m).^2 * Id + EIb * S)^(-1) ;
    y(m) = Wzc.' * temp * Wzc ;
end
end

function y = Psdef(om,T0, EI, rA, L)
lm      = sqrt((sqrt(T0^2 + 4*EI*rA*om.^2) - T0)/2/EI) ;
lp      = sqrt((T0 + sqrt(T0^2 + 4*EI*rA*om.^2))/2/EI) ;
y      = sin(lm*L) .* (1 + lm.^2 ./ lp.^2) ;
end

function y = Pspdef(om,T0, EI, rA, L)
lm      = sqrt((sqrt(T0^2 + 4*EI*rA*om.^2) - T0)/2/EI) ;
lp      = sqrt((T0 + sqrt(T0^2 + 4*EI*rA*om.^2))/2/EI) ;
y     = lm.*cos(lm*L) + (lm.^2.*sin(lm*L))./(lp.*tanh(lp*L)) ;
end

function y = Pspppdef(om,T0, EI, rA, L)
lm      = sqrt((sqrt(T0^2 + 4*EI*rA*om.^2) - T0)/2/EI) ;
lp      = sqrt((T0 + sqrt(T0^2 + 4*EI*rA*om.^2))/2/EI) ;
y   = -lm.^3.*cos(lm*L) + (lm.^2.*lp.*sin(lm*L))./(tanh(lp*L)) ;
end

function y = dPsdef(om,T0, EI, rA, L)
y       = sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)).*((2*EI*om*rA)./((T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2).*(T0^2 + 4*EI*rA*om^2).^(1/2)) + (2*EI*om*rA*(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2))./((T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2)^2*(T0^2 + 4*EI*rA*om^2).^(1/2))) - (L*om*rA*cos(L*(-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)).*((T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)./(T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2) - 1))./((-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2).*(T0^2 + 4*EI*rA*om^2).^(1/2)) ;
end

function y = dPspdef(om,T0, EI, rA, L)
y      = (om*rA*cos(L*(-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)))./((-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2).*(T0^2 + 4*EI*rA*om^2).^(1/2)) - (L*om*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)))./(T0^2 + 4*EI*rA*om^2).^(1/2) + (2*om*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)))./(tanh(L*((T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)).*((T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2).*(T0^2 + 4*EI*rA*om^2).^(1/2)) + (om*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)).*(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2))./(EI*tanh(L*((T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)).*((T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(3/2).*(T0^2 + 4*EI*rA*om^2).^(1/2)) - (L*om*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)).*(tanh(L*((T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2))^2 - 1).*(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2))./(tanh(L*((T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2))^2*(T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2).*(T0^2 + 4*EI*rA*om^2).^(1/2)) - (L*om*rA*cos(L*(-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)).*(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2))./(EI*tanh(L*((T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)).*(-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2).*((T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2).*(T0^2 + 4*EI*rA*om^2).^(1/2)) ;
end

function y = dPspppdef(om,T0, EI, rA, L)
y    = (2*om*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)).*((T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2))./(tanh(L*((T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)).*(T0^2 + 4*EI*rA*om^2).^(1/2)) - (3*om*rA*cos(L*(-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)).*(-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2))./(T0^2 + 4*EI*rA*om^2).^(1/2) - (L*om*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)).*(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2))./(EI*(T0^2 + 4*EI*rA*om^2).^(1/2)) - (om*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)).*(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2))./(EI*tanh(L*((T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)).*((T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2).*(T0^2 + 4*EI*rA*om^2).^(1/2)) - (L*om*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)).*(tanh(L*((T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2))^2 - 1).*(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2))./(EI*tanh(L*((T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2))^2*(T0^2 + 4*EI*rA*om^2).^(1/2)) - (L*om*rA*cos(L*(-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)).*(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2).*((T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2))./(EI*tanh(L*((T0/2 + (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2)).*(-(T0/2 - (T0^2 + 4*EI*rA*om^2).^(1/2)/2)/EI).^(1/2).*(T0^2 + 4*EI*rA*om^2).^(1/2)) ;
end

