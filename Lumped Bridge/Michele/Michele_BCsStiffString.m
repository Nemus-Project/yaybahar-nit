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
M     = 0.01 ;             %-- mass [Kg]
K     = 1e6 ;              %-- stiffness [N/m]

%--------------------------------------------------------

%-- derived paramters
A     = pi*r^2 ;
I     = pi*r^4/4 ;
EI    = E*I ;
rA    = rho*A ;



%--------------------------------------------------------

%-- newton-raphson routine for eigenvalues (eigenfrequencies)

om    = 20 * 2 * pi : 10000 * 2 * pi ;

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
omegaVec = zeros(length(zom),1) ;

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
    
    omegaVec(n) = om ;
    
    line([om,om],[-mm,mm],'linestyle','--','color','k') ;
    
end

%--------------------------------------------------------
%EigenModes
syms x lambdap lambdam cc
psi = (sin(lambdam*x) + cc*sinh(lambdap*x))^2;
Psi = int(psi,x);

rPlus    = (T0 + sqrt(T0^2 + 4*EI*rA*omegaVec.^2))/2/EI;
rMinus    = (sqrt(T0^2 + 4*EI*rA*omegaVec.^2) - T0)/2/EI;

lambdaM   = sqrt(rMinus);
lambdaP   = sqrt(rPlus);

%--------------------------------------------------------
%Checking modes orthogonality
lambdaP1 = lambdaP(10); lambdaM1 = lambdaM(10);
lambdaP2 = lambdaP(11); lambdaM2 = lambdaM(11);

syms y lambdam1 lambdam2 lambdap1 lambdap2 cc1 cc2 nc1 nc2
prod = nc1*(sin(lambdam1*y) + cc1*sinh(lambdap1*y)) * nc2*(sin(lambdam2*y) + cc2*sinh(lambdap2*y));
intProd = int(prod,y);
 
normCoeff1 = ComputeNormalizationFac(lambdaM1, lambdaP1,L);
normCoeff2 = ComputeNormalizationFac(lambdaM2, lambdaP2,L);
ip1 = InnerProd(L, lambdaM1, lambdaM2, lambdaP1, lambdaP2, L, normCoeff1, normCoeff2);
ip2 = InnerProd(0, lambdaM1, lambdaM2, lambdaP1, lambdaP2, L, normCoeff1, normCoeff2);
innerProduct = (ip1 - ip2);

 



%--------------------------------------------------------
%Functions
function mode = ComputeModes(x, lambdam, lambdap, L)
    sinhL    = sinh(lambdap*L);
    sinL    = sin(lambdam*L);

    const = (lambdam.^2./lambdap.^2).*(sinL./sinhL);
    mm = sin(lambdam*x) + const.*sinh(lambdap*x);

    D = ComputeNormalizationFac(lambdam, lambdap, L);
    mode = mm*D;
end

function norm = ComputeNormalizationFac(lambdam, lambdap, L)
    intL = IntegralPsiSq(L, lambdam, lambdap, L);
    intZero = IntegralPsiSq(0, lambdam, lambdap, L);

    norm = 1/sqrt(intL - intZero);
end

function int = IntegralPsiSq(x, lambdam, lambdap, L)
    sinhL    = sinh(lambdap*L);
    sinL    = sin(lambdam*L);

    const = (lambdam.^2./lambdap.^2).*(sinL./sinhL);
    int = x/2 - (const.^2*x)/2 - (cos(lambdam*x).*sin(lambdam*x))./(2*lambdam) + (const.^2.*cosh(lambdap*x).*sinh(lambdap*x))./(2*lambdap) - (2*const.*lambdam.*cos(lambdam*x).*sinh(lambdap*x))./(lambdam.^2 + lambdap.^2) + (2*const.*lambdap.*cosh(lambdap*x).*sin(lambdam*x))./(lambdam.^2 + lambdap.^2);
end

function int = InnerProd(x, lambdam1, lambdam2, lambdap1, lambdap2, L, nc1, nc2)
    cc1 = (lambdam1^2/lambdap1^2)*(sin(lambdam1*L)/sinh(lambdap1*L));
    cc2 = (lambdam2^2/lambdap2^2)*(sin(lambdam2*L)/sinh(lambdap2*L));
     
    int = (lambdam2*nc1*nc2*cos(lambdam2*x)*sin(lambdam1*x))/(lambdam1^2 - lambdam2^2) - (lambdam1*nc1*nc2*cos(lambdam1*x)*sin(lambdam2*x))/(lambdam1^2 - lambdam2^2) - (cc1*lambdam2*nc1*nc2*cos(lambdam2*x)*sinh(lambdap1*x))/(lambdam2^2 + lambdap1^2) - (cc2*lambdam1*nc1*nc2*cos(lambdam1*x)*sinh(lambdap2*x))/(lambdam1^2 + lambdap2^2) + (cc1*lambdap1*nc1*nc2*cosh(lambdap1*x)*sin(lambdam2*x))/(lambdam2^2 + lambdap1^2) + (cc2*lambdap2*nc1*nc2*cosh(lambdap2*x)*sin(lambdam1*x))/(lambdam1^2 + lambdap2^2) + (cc1*cc2*lambdap1*nc1*nc2*cosh(lambdap1*x)*sinh(lambdap2*x))/(lambdap1^2 - lambdap2^2) - (cc1*cc2*lambdap2*nc1*nc2*cosh(lambdap2*x)*sinh(lambdap1*x))/(lambdap1^2 - lambdap2^2);
end
        
        
        
        
        
        
        
        
        
