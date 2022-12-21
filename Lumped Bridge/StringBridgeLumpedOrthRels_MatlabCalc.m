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

rPlus    = (T0 + sqrt(T0^2 + 4*EI*rA*omegaVec.^2))/2/EI;
rMinus    = (sqrt(T0^2 + 4*EI*rA*omegaVec.^2) - T0)/2/EI;

lambdaM   = sqrt(rMinus);
lambdaP   = sqrt(rPlus);

%inner product with different eigenmodes (should be ~0)
inProd1 = CheckOrthogonality(lambdaM(30), lambdaM(40), lambdaP(30), lambdaP(40), L)
%inner product with same eigenmodes (should be 1)
inProd2 = CheckOrthogonality(lambdaM(20), lambdaM(20), lambdaP(20), lambdaP(20), L)
%critical modes
inProd3 = CheckOrthogonality(lambdaM(11), lambdaM(12), lambdaP(11), lambdaP(12), L)

%--------------------------------------------------------
function value = ComputeNormalizationFac(lambdam, lambdap, L)
    sinhL    = sinh(lambdap*L);
    sinL    = sin(lambdam*L);
    const = (lambdam.^2./lambdap.^2).*(sinL./sinhL);

    syms x
    psi = (sin(lambdam*x) + const.*sinh(lambdap*x)).*(sin(lambdam*x) + const.*sinh(lambdap*x));
    integral = int(psi, 0, L);
    Fvpa = vpa(integral);
    value = 1/sqrt(Fvpa);
end

function value = CheckOrthogonality(lambdam1, lambdam2, lambdap1, lambdap2, L)
    sinhL1    = sinh(lambdap1*L);
    sinL1    = sin(lambdam1*L);
    const1 = (lambdam1.^2./lambdap1.^2).*(sinL1./sinhL1);

    sinhL2    = sinh(lambdap2*L);
    sinL2    = sin(lambdam2*L);
    const2 = (lambdam2.^2./lambdap2.^2).*(sinL2./sinhL2);

    normCoeff1 = double(ComputeNormalizationFac(lambdam1, lambdap1, L));
    normCoeff2 = double(ComputeNormalizationFac(lambdam2, lambdap2, L));

    syms x
    psi = normCoeff1*normCoeff2*(sin(lambdam1*x) + const1*sinh(lambdap1*x))*(sin(lambdam2*x) + const2*sinh(lambdap2*x));

    integral = vpaintegral(psi, 0, L);
    value = integral;
end

