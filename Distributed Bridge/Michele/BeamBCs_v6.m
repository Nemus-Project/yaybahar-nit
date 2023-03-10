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

% r = 3.75e-04;                  % string radius [m]
%         rho = 3.7575e3;                     % string Density [kg/m] 
%         T0 = 153;                           % tension [N] 
%         E = 25e9;
%         L = 0.69;



%-- bridge parameters
Eb    = 3e10 ;             %-- Young's mod [Pa]
rhob  = 8e3 ;              %-- density [N/m]
rb    = 1e-3 ;             %-- radius [m]
Lb    = 7e-2 ;             %-- length [m]
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

% while m < 1
%     m = m + 1 ;
%     Om = sqrt(EIb/rAb) * (m^2 * pi^2 / Lb^2) ;
%     Wzc = [Wzc; sqrt(2/Lb)*sin(m*pi*zc)] ;
%     S   = [S; m^4*pi^4/Lb^4] ;
% end

Mmax   = m ;

S = diag(S) ;
Id = eye(Mmax) ;

% Om = Om(1) ;
% Wzc = Wzc(1) ;
% S = S(1) ;



%--------------------------------------------------------

%-- newton-raphson routine for eigenvalues

om    = 1 * 2 * pi : MaxOm ;
g     = gdef(om,rAb,EIb,Wzc,S) ;
gp    = gpdef(om,rAb,EIb,Wzc,S) ;
Ps    = Psdef(om,T0, EI, rA, L) ;
Psp   = Pspdef(om,T0, EI, rA, L) ;
Psppp = Pspppdef(om,T0, EI, rA, L) ;

% %-- eigenvale equation
fom   =  Ps - g.* (-T0*Psp + EI*Psppp) ;
plot(om,fom) ; hold on ; grid on ;
mm = max(abs(fom)) ;
ylim([-10,10])



%-- search for zeros (approximate search)
zom = [] ;
for n = 1 : length(om) - 1

    tt = fom(n) * fom(n+1) ;

    if tt <= 0

        if abs(fom(n)) < 1

            zom = [zom;om(n)] ;

        end

    end

end


%-- initialise zero vector
zerF = zeros(length(zom),1) ;

%-- newton Raphson
for n = 1 : length(zom)

    om = zom(n) ;
    tol = 1 ;
    iter = 0 ;

    while tol > 1e-14 && iter < 200

        iter  = iter + 1 ;

        g        = gdef(om,rAb,EIb,Wzc,S) ;
        gp       = gpdef(om,rAb,EIb,Wzc,S) ;
        Ps       = Psdef(om,T0, EI, rA, L) ;
        Psp      = Pspdef(om,T0, EI, rA, L) ;
        Psppp    = Pspppdef(om,T0, EI, rA, L) ;
        dPs      = dPsdef(om,T0, EI, rA, L) ;
        dPsp     = dPspdef(om,T0, EI, rA, L) ;
        dPsppp   = dPspppdef(om,T0, EI, rA, L) ;


        %         fom     = Ps - g.*(-T0*Psp + EI*Psppp) ;
        %         dfom    = dPs - gp.*(-T0*Psp + EI*Psppp) - g.*(-T0*dPsp + EI*dPsppp) ;
        fom     = 1./g .* Ps - (-T0*Psp + EI*Psppp) ;
        dfom    = -(gp./g.^2) .* Ps + 1./g .* dPs - (-T0*dPsp + EI*dPsppp) ;
        omk     = om - fom/dfom ;

        tol     = abs(omk-om) ;
        om      = omk ;

    end

    %     g        = gdef(om,rAb,EIb,Wzc,S) ;
    %     Ps       = Psdef(om,T0, EI, rA, L) ;
    %     Psp      = Pspdef(om,T0, EI, rA, L) ;
    %     Psppp    = Pspppdef(om,T0, EI, rA, L) ;
    %     fom      = Ps - g.*(-T0*Psp + EI*Psppp) ;
    %     zerF(n) = [0] ;
    %     if abs(fom) < 1
    zerF(n) = om ;
    line([om,om],[-mm,mm],'linestyle','--','color','k') ;

    %     end

end

omSS = zeros(length(zom),1) ;
for m = 1 : length(zom)

    omSS(m) = sqrt(T0/rA * (m*pi/L)^2 + EI/rA * (m*pi/L)^4);
end


figure;

plot(omSS,'*') ; hold on; plot(zerF,'+') ;



%---- build mass and stiffness matrices
Nmodes      = 70 ;
bvec        = zeros(Nmodes,1) ;
for n = 1 : Nmodes
    temp = Wzc.' * (zerF(n)^2 * rAb * Id - EIb * S)^(-1) * Wzc ;
    bvec(n) = psiDef(T0,L,sqrt(EI),rA,zerF(n),L) / temp ;
end


% dx          = 1e-5 * L ;
% x           = 0 : dx : L ;

% if mod(length(x),2) == 0
%     x = [x,x(end)+dx] ;
%     x = x*L/(L+dx) ;
%     dx = x(2) - x(1) ;
% end
%
% simpWeights = 4 * ones(1,length(x)) ;
% for n = 2 : length(x) - 2
%     if mod(n,2) == 0
%         simpWeights(n+1) = 2 ;
%     end
% end
% simpWeights(1)      = 1 ;
% simpWeights(2)      = 4 ;
% simpWeights(end)    = 1 ;
% simpWeights         = simpWeights/3 ;
%
% omvec       = zerF(1:Nmodes) ;
%
% psiMat      = zeros(Nmodes,length(x)) ;
% psipMat     = zeros(Nmodes,length(x)) ;
% psippMat    = zeros(Nmodes,length(x)) ;
% psipppMat   = zeros(Nmodes,length(x)) ;
%
%
%
% for n = 1 : Nmodes
%
%     om              = omvec(n) ;
%     kappa           = sqrt(EI) ;
%     psiMat(n,:)     = psiDef(T0,L,kappa,rA,om,x) ;
%     psipMat(n,:)    = psipDef(T0,L,kappa,rA,om,x) ;
%     psippMat(n,:)   = psippDef(T0,L,kappa,rA,om,x) ;
%     psipppMat(n,:)  = psipppDef(T0,L,kappa,rA,om,x) ;
%
% end


MassMat      = zeros(Nmodes,Nmodes) ;
StiffMat     = zeros(Nmodes,Nmodes) ;

% gg           = (S * Wzc) ;
% fInv         = gg * gg.' ;
wMass        = rAb / (Wzc.' * Wzc) ;
wStiff       = EIb / (Wzc.' * Wzc) ;
PsiL         = zeros(Nmodes,1) ;

for n = 1 : Nmodes
   PsiL(n) =  psiDef(T0,L,sqrt(EI),rA,zerF(n),L) ;
end

for m = 1 : Nmodes
    for n = m : Nmodes
        lmm      = sqrt((sqrt(T0^2 + 4*EI*rA*zerF(m)^2) - T0)/2/EI) ;
        lpm      = sqrt((T0 + sqrt(T0^2 + 4*EI*rA*zerF(m)^2))/2/EI) ;
        lmn      = sqrt((sqrt(T0^2 + 4*EI*rA*zerF(n)^2) - T0)/2/EI) ;
        lpn      = sqrt((T0 + sqrt(T0^2 + 4*EI*rA*zerF(n)^2))/2/EI) ;
        MassMat(m,n)   = rA * MassMatInt(m,n,lmm,lmn,lpm,lpn,L)  ;
        % MassMat2(m,n)  = rA * sum(psiMat(m,:).*psiMat(n,:).*simpWeights) * dx  ; %+ wMass*psiMat(m,end).*psiMat(n,end) ;
        % StiffMat(m,n)  = T0 * sum(psipMat(m,:).*psipMat(n,:).*simpWeights) * dx + EI * sum(psippMat(m,:).*psippMat(n,:).*simpWeights) * dx ; % + wStiff*psiMat(m,end).*psiMat(n,end) ;
        StiffMat(m,n)  = T0  * StiffMatIntPrime(m,n,lmm,lmn,lpm,lpn,L) + EI * StiffMatIntPrimePrime(m,n,lmm,lmn,lpm,lpn,L)  ;
        if m == n
            MassMat(m,n)  = 0.5 * MassMat(m,n) ;
            StiffMat(m,n) = 0.5 * StiffMat(m,n) ;
        end
    end
end

MassMat  = MassMat + MassMat.' ;
%MassMat  = diag(diag(MassMat)) ;
StiffMat = StiffMat + StiffMat.' ;

for m = 1 : Nmodes
    for n = 1 : Nmodes
        %  StiffMat(m,n)  = StiffMat(m,n) - psiDef(T0,L,sqrt(EI),rA,zerF(m),L) * (T0 * psipDef(T0,L,sqrt(EI),rA,zerF(n),L) - EI * psipppDef(T0,L,sqrt(EI),rA,zerF(n),L))  ;
        StiffMat(m,n)  = StiffMat(m,n) - PsiL(m) * bvec(n) ;
    end
end


% SMat = zeros(Nmodes+Mmax,Nmodes+Mmax) ;
% MMat = zeros(Nmodes+Mmax,Nmodes+Mmax) ;
% 
% SMat(1:Nmodes,1:Nmodes) = StiffMat - PsiL * bvec.' ;
% MMat(1:Nmodes,1:Nmodes) = MassMat ;
% MMat(end-Mmax+1:end,end-Mmax+1:end) = rAb * Id ;
% SMat(1:Nmodes,end-Mmax+1:end) = 0 * wStiff * PsiL * (Wzc.' * S) ;
% SMat(end-Mmax+1:end,end-Mmax+1:end) = EIb * S ;
% SMat(end-Mmax+1:end,1:Nmodes) = Wzc * bvec.' ;


[Vm,Dm] = eig(MassMat) ;
[V,D] = eig(MassMat \ StiffMat) ;
% [V,D] = eig(MMat \ SMat) ;


omm = sqrt(diag(D) ) ;
omm = sort(omm(1:Nmodes)) ;
err_percent = (omm(1:Nmodes)-zerF(1:Nmodes))./zerF(1:Nmodes)*100

eigenVals = eig(MassMat\StiffMat)

%-- project initial conditions

% u0 = sin(4*pi*x/L) ;
% projvec = zeros(Nmodes,1) ;
% for n = 1 : Nmodes
%     projvec(n) = psiMat(n,:) * u0.' * dx ;
% end
% q0 = MassMat \ (projvec) ;
% Kd = MassMat \ StiffMat ;

return
%--------------------------------------------------------------------------
function y = psiDef(T0,L,kappa,rA,om,x)

lp      = sqrt( (T0 + sqrt(T0^2 + 4*kappa^2*rA*om^2) ) / 2 / kappa^2 ) ;
lm      = sqrt( (sqrt(T0^2 + 4*kappa^2*rA*om^2) - T0 ) / 2 / kappa^2 ) ;
y       = sin(lm.*x) + lm.^2 .* sin(lm.*L) ./ lp.^2 ./ sinh(lp.*L) .* sinh(lp.*x) ;

end

function y = psipDef(T0,L,kappa,rA,om,x)

lp      = sqrt( (T0 + sqrt(T0^2 + 4*kappa^2*rA*om^2) ) / 2 / kappa^2 ) ;
lm      = sqrt( (sqrt(T0^2 + 4*kappa^2*rA*om^2) - T0 ) / 2 / kappa^2 ) ;
y       = lm*cos(lm.*x) + lm.^2 .* sin(lm.*L) ./ lp.^2 ./ sinh(lp.*L) * lp * cosh(lp.*x) ;

end

function y = psippDef(T0,L,kappa,rA,om,x)

lp      = sqrt( (T0 + sqrt(T0^2 + 4*kappa^2*rA*om^2) ) / 2 / kappa^2 ) ;
lm      = sqrt( (sqrt(T0^2 + 4*kappa^2*rA*om^2) - T0 ) / 2 / kappa^2 ) ;
y       = -lm^2*sin(lm.*x) + lm.^2 .* sin(lm.*L) ./ lp.^2 ./ sinh(lp.*L) * lp^2 * sinh(lp.*x) ;

end

function y = psipppDef(T0,L,kappa,rA,om,x)

lp      = sqrt( (T0 + sqrt(T0^2 + 4*kappa^2*rA*om^2) ) / 2 / kappa^2 ) ;
lm      = sqrt( (sqrt(T0^2 + 4*kappa^2*rA*om^2) - T0 ) / 2 / kappa^2 ) ;
y       = -lm^3*cos(lm.*x) + lm.^2 .* sin(lm.*L) ./ lp.^2 ./ sinh(lp.*L) * lp^3 * cosh(lp.*x) ;

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

function y = MassMatInt(m,n,lmm,lmn,lpm,lpn,L)

cm = lmm^2*sin(lmm*L)/lpm^2/sinh(lpm*L) ;
cn = lmn^2*sin(lmn*L)/lpn^2/sinh(lpn*L) ;

if m ~= n

    i1 = -(lmm*cos(L*lmm)*sin(L*lmn) - lmn*cos(L*lmn)*sin(L*lmm))/(lmm^2 - lmn^2) ;
    i2 = -(lmm*cos(L*lmm)*sinh(L*lpn) - lpn*cosh(L*lpn)*sin(L*lmm))/(lmm^2 + lpn^2) ;
    i3 = (lpm*cosh(L*lpm)*sinh(L*lpn) - lpn*cosh(L*lpn)*sinh(L*lpm))/(lpm^2 - lpn^2) ;
    i4 = -(lmn*cos(L*lmn)*sinh(L*lpm) - lpm*cosh(L*lpm)*sin(L*lmn))/(lmm^2 + lpn^2) ;

    y = i1 + cn*i2 + cm*i4 + cm*cn*i3 ;

else

    i4 = L/2 - sin(2*L*lmm)/(4*lmm) ;
    i5 = sinh(2*L*lpm)/(4*lpm) - L/2 ;
    i6 = -(cos(L*lpm)*sinh(L*lpm) - cosh(L*lpm)*sin(L*lpm))/(2*lpm) ;

    y = i4 + 2*cm*i6 + cm^2*i5 ;

end

end

function y = StiffMatIntPrime(m,n,lmm,lmn,lpm,lpn,L)

cm = lmm^2*sin(lmm*L)/lpm^2/sinh(lpm*L) ;
cn = lmn^2*sin(lmn*L)/lpn^2/sinh(lpn*L) ;

if m ~= n

    i1 = lmm*lmn*(sin(L*(lmm - lmn))/(2*(lmm - lmn)) + sin(L*(lmm + lmn))/(2*(lmm + lmn)))  ;
    i2 = lmm*lpn*cn*((lmm*cosh(L*lpn)*sin(L*lmm) + lpn*cos(L*lmm)*sinh(L*lpn))/(lmm^2 + lpn^2)) ;
    i3 = cm*cn*lpm*lpn*(lpm*cosh(L*lpn)*sinh(L*lpm) - lpn*cosh(L*lpm)*sinh(L*lpn))/(lpm^2 - lpn^2) ;
    i4 = lmn*lpm*cm*((lmn*cosh(L*lpm)*sin(L*lmn) + lpm*cos(L*lmn)*sinh(L*lpm))/(lmn^2 + lpm^2)) ;

    y = i1 + i2 + i3 + i4 ;

else

    i5 = lmm^2 * (L/2 + sin(2*L*lmm)/(4*lmm)) ;
    i6 = lpm^2 * cm^2 * (L/2 + sinh(2*L*lpm)/(4*lpm)) ;
    i7 = 2 * cm * lpm * lmm * (sin(L*lpm*(1 - 1i))*(1 + 1i) + sin(L*lpm*(1 + 1i))*(1 - 1i))/(4*lpm) ;

    y = i5 + i6 + i7 ;

end

end


function y = StiffMatIntPrimePrime(m,n,lmm,lmn,lpm,lpn,L)

cm = lmm^2*sin(lmm*L)/lpm^2/sinh(lpm*L) ;
cn = lmn^2*sin(lmn*L)/lpn^2/sinh(lpn*L) ;

if m ~= n

    i1 = -(lmm*cos(L*lmm)*sin(L*lmn) - lmn*cos(L*lmn)*sin(L*lmm))/(lmm^2 - lmn^2) ;
    i2 = -(lmm*cos(L*lmm)*sinh(L*lpn) - lpn*cosh(L*lpn)*sin(L*lmm))/(lmm^2 + lpn^2) ;
    i3 = (lpm*cosh(L*lpm)*sinh(L*lpn) - lpn*cosh(L*lpn)*sinh(L*lpm))/(lpm^2 - lpn^2) ;
    i4 = -(lmn*cos(L*lmn)*sinh(L*lpm) - lpm*cosh(L*lpm)*sin(L*lmn))/(lmm^2 + lpn^2) ;

    y = lmm^2 * lmn^2 * i1 - lmm^2 * lpn^2 * cn * i2 - lmn^2 * lpm^2 * cm * i4 + lpm * lpn * cm * cn * i3 ;

else

    i5 = L/2 - sin(2*L*lmm)/(4*lmm) ;
    i6 = sinh(2*L*lpm)/(4*lpm) - L/2 ;
    i7 = -(cos(L*lpm)*sinh(L*lpm) - cosh(L*lpm)*sin(L*lpm))/(2*lpm) ;

    y = lmm^4 * i5 - 2*cm*lpm^2*lmm^2*i7 + lpm^4*cm^2*i6 ;

end

end

