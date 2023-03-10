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
% T0    = 147.7 ;            %-- tension [N]
% r     = 4.55e-4 ;          %-- radius [m]
% E     = 2.5e10 ;           %-- Young's mod [Pa]
% rho   = 5535 ;             %-- density [kg/m^3]
L     = 0.69 ;             %-- length [m]
r = 3.75e-04;                  % string radius [m]
        rho = 3.7575e3;                     % string Density [kg/m] 
        T0 = 153;                           % tension [N] 
        E = 25e9;                          % young modulus [Pa]


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
while Om<=MaxFreq*2*pi%m <= 10
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

% %-- eigenvale equation
 fom   =  Ps - g.* (-T0*Psp + EI*Psppp) ;
 plot(om,fom) ; hold on ; grid on ;
 mm = max(abs(fom)) ;
 ylim([-10,10])



%-- search for zeros (approximate search)
zom = [] ;
zomNoCheck = [];
for n = 1 : length(om) - 1

    tt = fom(n) * fom(n+1) ;

    if tt <= 0

        if abs(fom(n)) < 1

        zom = [zom;om(n)] ;

        end

        zomNoCheck = [zomNoCheck;om(n)] ;

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
%           fom     = 1./g .* Ps - (-T0*Psp + EI*Psppp) ;
%           dfom    = -(gp./g.^2) .* Ps + 1./g .* dPs - (-T0*dPsp + EI*dPsppp) ;
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

figure;

plot(zom,'*') ; hold on; plot(zerF,'+') ;


%Confronto tra eigenfreqs simply supported e non
% omSS = zeros(length(zom),1) ;
% for m = 1 : length(zom)
% 
%     omSS(m) = sqrt(T0/rA * (m*pi/L)^2 + EI/rA * (m*pi/L)^4);
% end


% figure;
% 
% plot(omSS,'*') ; hold on; plot(zerF,'+') ;



%---- build mass and stiffness matrices
% Nmodes      = 10 ;
Nmodes = length(zerF);
dx          = 1e-5 * L ;
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


MassMat      = zeros(Nmodes,Nmodes) ;
StiffMat     = zeros(Nmodes,Nmodes) ;

for m = 1 : Nmodes

    for n = 1 : Nmodes 

        MassMat(m,n)  = rA * (psiMat(m,1:end-1)*(psiMat(n,1:end-1).') + 0.5 * psiMat(m,end)*(psiMat(n,end))) * dx ;
        StiffMat(m,n) = T0 * (psipMat(m,1:end-1)*(psipMat(n,1:end-1).') + 0.5 * psipMat(m,end)*(psipMat(n,end)) ) * dx + EI * (psippMat(m,:)*(psippMat(n,:).') + 0.5 * psippMat(m,end)*(psippMat(n,end)) ) * dx - psiMat(m,end) * (T0 * psipMat(n,end) - EI * psipppMat(n,end)) ;
    
    end

end

[V,D] = eig(StiffMat,MassMat) ;
% disp(['MassMat symmetry', string(issymmetric(MassMat))]);
% disp(['StiffMat symmetry', string(issymmetric(StiffMat))]);

omm = sort(sqrt(diag(D))) ;

%checks if the eigenfreqs calculated with newton raphson and the ones
%obtained as eigenvalues of mass and stiff matrices are almost equal
err_percent = (omm(1:Nmodes)-zerF(1:Nmodes))./zerF(1:Nmodes)*100

%% Simulation
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Initializing oscillators and coefficients
SR = 44100;
k = 1/SR;
timeSamples = SR;

for i = 1:Nmodes
    if k>=(2/omm(i))
        disp("stability condition violated");
        %return;
    end
end

% q = zeros(modesNumber,1);
q = ones(Nmodes,1);
qPrev = q;
qNext = zeros(Nmodes,1);

qDense = q; qDensePrev = qPrev; qDenseNext = qNext;

sigma0 = 0;

%Qua c'è un confronto per vedere cosa cambia se simulo il sistema con le
%autofrequenze (autovalori di MassMat^-1*StiffMat o zeri della charact
%equation, ho visto che questi 2 corrispondono con la verifica in
%err_percent) o direttamente con l'inversione di matrice. Ho visto che la
%matrice MassMat^-1*StiffMat non è diagonale, ma quasi, nel senso che gli
%elementi sulla diagonale sono 2 o più odg più grandi degli altri, quando
%poi si va a simulare moltiplicando per k^2 gli elementi fuori dalla
%diagonale praticamente si annullano. Ma soprattutto, gli elementi sulla
%diagonale sono esattamente gli autivalori! Questa cosa è strana, quindi in
%pratica la diagonale di MassMat^-1*StiffMat e i suoi autovalori
%coincidono! Quindi non ho bisogno di un cambio di variabile alla fine per
%diagonalizzare il sistema. err_percent 2 fa vedere la differenza tra la
%diagonale di MassMat^-1*StiffMat e i suoi autvalori. Il punto è che dato
%che MassMat non è simmetrica, non è detto che abbia autovalori reali!
invMassMat = MassMat^-1;
RsMat = invMassMat*StiffMat;
err_percent2 = (sort(diag(RsMat)) - sort(diag(D)))./sort(diag(D)) * 100

Mat = 2 - omm.^2*k^2;
MatDense = 2 - RsMat*k^2;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Simulation
close all
for n = 1:timeSamples
    hold off
    qNext = q.*Mat - qPrev;

    qPrev = q;
    q = qNext;
    
    %devo prendere solo la diagonale però, altrimenti se faccio
    %direttamente l'inversione di matrice esplode, non so perchè
%     qDenseNext = MatDense*qDense - qDensePrev;
    qDenseNext = diag(MatDense).*qDense - qDensePrev;
 
    qDensePrev = qDense;
    qDense = qDenseNext;

    plot(q,'*');
    hold on
    plot(qDense);
    ylim([-20,20]);
    drawnow
end

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


