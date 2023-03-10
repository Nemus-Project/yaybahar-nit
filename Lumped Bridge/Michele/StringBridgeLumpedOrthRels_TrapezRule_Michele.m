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
M     = 1e-2 ;             %-- mass [Kg]
K     = 1e6 ;             %-- stiffness [N/m]

%--------------------------------------------------------

%-- derived paramters
A     = pi*r^2 ;
I     = pi*r^4/4 ;
EI    = E*I ;
rA    = rho*A ;



%--------------------------------------------------------

%-- newton-raphson routine for eigenvalues (eigenfrequencies)

om    = 20 * 2 * pi : 0.1 : 10000 * 2 * pi ;

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
figure(1)
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
omegaVec = zeros(length(zom),1) ;

%-- newton Raphson
for n = 1 : length(zom)
    
    om = zom(n) ;
    tol = 1 ;
    iter = 0 ;
    
    while tol > 1e-13 && iter < 20
        
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
    
    omegaVec(n) = om ;
    
    line([om,om],[-mm,mm],'linestyle','--','color','k') ;
    
end



%--------------------------------------------------------
%EigenModes

dx          = 1e-6 * L ;
x           = (0 : dx : L).' ;
Nmod        = 20 ;
massMatrix  = zeros(Nmod,Nmod) ;
stiffMatrix  = zeros(Nmod,Nmod) ;


for m = 1 : Nmod
    for n = m : Nmod
        
        omm = omegaVec(m) ;
        omn = omegaVec(n) ;
        
        rPm   = (T0 + sqrt(T0^2 + 4*EI*rA*omm^2))/2/EI;
        rMm   = (sqrt(T0^2 + 4*EI*rA*omm^2) - T0)/2/EI;
        lMm   = sqrt(rMm);
        lPm   = sqrt(rPm);
        
        rPn   = (T0 + sqrt(T0^2 + 4*EI*rA*omn^2))/2/EI;
        rMn   = (sqrt(T0^2 + 4*EI*rA*omn^2) - T0)/2/EI;
        lMn   = sqrt(rMn);
        lPn   = sqrt(rPn);
        
        shLm  = sinh(lPm*L);
        snLm  = sin(lMm*L);
        cm    = (lMm^2/lPm^2)*(snLm/shLm);
        
        shLn  = sinh(lPn*L);
        snLn  = sin(lMn*L);
        cn    = (lMn^2/lPn^2)*(snLn/shLn);
        
        fm    = sin(lMm*x) + cm*sinh(lPm*x) ;
        fn    = sin(lMn*x) + cn*sinh(lPn*x) ;

        fpm = lMm*cos(lMm*x) + cm*lPm*cosh(lPm*x);
        fpn = lMn*cos(lMn*x) + cn*lPn*cosh(lPn*x);

        fppm = -lMm^2*sin(lMm*x) + cm*lPm^2*sinh(lPm*x);
        fppn = -lMn^2*sin(lMn*x) + cn*lPn^2*sinh(lPn*x);

        fmL   = sin(lMm*L) + cm*sinh(lPm*L) ;
        fnL   = sin(lMn*L) + cn*sinh(lPn*L) ;
        
        %-- trapezoid rule
        fm(1) = 0.5*fm(1) ; fm(end) = 0.5*fm(end) ;
        fpm(1) = 0.5*fpm(1) ; fpm(end) = 0.5*fpm(end) ;
        fppm(1) = 0.5*fppm(1) ; fppm(end) = 0.5*fppm(end) ;

        massMatrix(m,n) = rA*fm.'*fn*dx + M*fmL*fnL;
        stiffMatrix(m,n) = (T0*(fpm.'*fpn) + EI*(fppm.'*fppn))*dx + K*fmL*fnL;

        %l'ordine delle autofunzioni non cambia quindi le matrici sono
        %simmetriche
        if m ~= n
            massMatrix(n,m) = massMatrix(m,n) ;
            stiffMatrix(n,m) = stiffMatrix(m,n);
        end
        
    end
end

figure(2)
imagesc(massMatrix);
figure(3)
imagesc(stiffMatrix);

omegaSqCalc = massMatrix^-1*stiffMatrix;

omegaDiff = (diag(omegaSqCalc) - omegaVec(1:Nmod).^2)./omegaVec(1:Nmod).^2;

