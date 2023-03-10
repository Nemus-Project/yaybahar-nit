clear all
close all

%--------------------------------------------------------
%          Modes of a String Terminated with
%               Mass - Spring
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
M     = 0.07 ;             %-- mass [Kg]
K     = 1e4 ;              %-- stiffness [N/m]

%--------------------------------------------------------

%-- derived paramters
A     = pi*r^2 ;
I     = pi*r^4/4 ;
EI    = E*I ;
rA    = rho*A ;



%--------------------------------------------------------

%-- newton-raphson routine for eigenvalues

om    = 1 * 2 * pi : 5000 * 2 * pi ;

rp    = (T0 + sqrt(T0^2 + 4*EI*rA*om.^2))/2/EI ;
rm    = (sqrt(T0^2 + 4*EI*rA*om.^2) - T0)/2/EI ;

gm1   = sqrt(rm) ;
gm2   = sqrt(rp) ;

th    = tanh(gm2*L) ;
cs    = cos(gm1*L) ;
sn    = sin(gm1*L) ;
g     = M*(-om.^2) + K  ;

%-- eigenvale equation
fom   = -T0*(gm1.*th.*cs + gm1.^2./gm2.*sn) + EI*(-gm1.^3.*th.*cs + gm1.^2.*gm2.*sn) - g.*th.*sn.*(1+gm1.^2./gm2.^2) ;
plot(om,fom) ; hold on ; grid on ;
mm = max(abs(fom)) ;


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
        rp    = (T0 + sqrt(T0^2 + 4*EI*rA*om.^2))/2/EI ;
        rm    = (sqrt(T0^2 + 4*EI*rA*om.^2) - T0)/2/EI ;
        gm1   = sqrt(rm) ;
        gm2   = sqrt(rp) ;
        th    = tanh(gm2*L) ;
        cs    = cos(gm1*L) ;
        sn    = sin(gm1*L) ;
        g     = M*(-om.^2) + K  ;

        fom     = -T0*(gm1.*th.*cs + gm1.^2./gm2.*sn) + EI*(-gm1.^3.*th.*cs + gm1.^2.*gm2.*sn) - g.*th.*sn.*(1+gm1.^2./gm2.^2) ;
        dfom    = (L.*om.*rA.*cos(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*(- M.*om.^2 + K).*((T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./(T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2) - 1))./((-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2).*(T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)) - T0.*((2.*om.*rA.*sin(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)))./(((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2).*(T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)) + (om.*rA.*cos(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)))./((-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2).*(T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)) - (L.*om.*rA.*sin(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)))./(T0.^2 + 4.*EI.*rA.*om.^2).^(1./2) + (om.*rA.*sin(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2))./(EI.*((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(3./2).*(T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)) - (L.*om.*rA.*cos(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*(tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).^2 - 1).*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2))./(((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2).*(T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)) - (L.*om.*rA.*cos(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2))./(EI.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2).*((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2).*(T0.^2 + 4.*EI.*rA.*om.^2).^(1./2))) - sin(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*(- M.*om.^2 + K).*((2.*EI.*om.*rA)./((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2).*(T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)) + (2.*EI.*om.*rA.*(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2))./((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2).^2.*(T0.^2 + 4.*EI.*rA.*om.^2).^(1./2))) - 2.*M.*om.*sin(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*((T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./(T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2) - 1) - (L.*om.*rA.*sin(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*(tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).^2 - 1).*(- M.*om.^2 + K).*((T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./(T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2) - 1))./(((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2).*(T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)) - EI.*((3.*om.*rA.*cos(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2))./(T0.^2 + 4.*EI.*rA.*om.^2).^(1./2) - (2.*om.*rA.*sin(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2))./(T0.^2 + 4.*EI.*rA.*om.^2).^(1./2) + (om.*rA.*sin(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2))./(EI.*((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2).*(T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)) - (L.*om.*rA.*cos(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*(tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).^2 - 1).*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(3./2))./(((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2).*(T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)) + (L.*om.*rA.*sin(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2))./(EI.*(T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)) + (L.*om.*rA.*cos(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2)).*(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2).*((T0./2 + (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2))./(EI.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*om.^2).^(1./2)./2)./EI).^(1./2).*(T0.^2 + 4.*EI.*rA.*om.^2).^(1./2))) ;

        omk     = om - fom/dfom ;

        tol     = abs(omk-om) ;
        om      = omk ;

    end

    zerF(n) = om ;

    line([om,om],[-mm,mm],'linestyle','--','color','k') ;

end

zerF = sort(zerF) ;


%---- build mass and stiffness matrices
Nmodes      = 10 ;
dx          = 1e-6 * L ;
x           = 0 : dx : L ;
omvec       = zerF(1:Nmodes) ;

psiMat      = zeros(Nmodes,length(x)) ;
psipMat     = zeros(Nmodes,length(x)) ;
psippMat    = zeros(Nmodes,length(x)) ;
psipppMat   = zeros(Nmodes,length(x)) ;



for n = 1 : Nmodes

    om              = omvec(n) ;
    kappa           = sqrt(EI) ;
    psiMat(n,:)     = psiDef(T0,L,kappa,rA,om,x) ;
    psipMat(n,:)    = psipDef(T0,L,kappa,rA,om,x) ;
    psippMat(n,:)   = psippDef(T0,L,kappa,rA,om,x) ;
    psipppMat(n,:)  = psipppDef(T0,L,kappa,rA,om,x) ;

end


MassMat = zeros(Nmodes,Nmodes) ;
StiffMat = zeros(Nmodes,Nmodes) ;

for m = 1 : Nmodes

    for n = 1 : Nmodes 

        MassMat(m,n) = rA * (psiMat(m,:)*(psiMat(n,:).')) * dx ;
        StiffMat(m,n) = T0 * (psipMat(m,:)*(psipMat(n,:).')) * dx + EI * (psippMat(m,:)*(psippMat(n,:).')) * dx - psiMat(m,end) * (T0 * psipMat(n,end) - EI * psipppMat(n,end)) ;
    
    end

end
det(MassMat)
[V,D] = eig(StiffMat,MassMat) ;
inv(MassMat);

omm = sort(sqrt(diag(D))) ;


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








