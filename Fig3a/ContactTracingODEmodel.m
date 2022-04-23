function udot= ContactTracingODEmodel(t,u,par)

uS=u(1);
uE=u(2);
uA=u(3);
uI=u(4);
uR=u(5);
uTs=u(6);
uTr=u(7);
uTq=u(8);
uQ=u(9);
%uRq=u(10);

k=par.k;
beta=par.beta;
N=par.N;
tau=par.tau;
delta=par.delta;
alphaA=par.alphaA;
alphaI=par.alphaI;
gammaA=par.gammaA;
gammaI=par.gammaI;
sigma=par.sigma;
gammaQ=par.gammaQ;

bhat=k*beta/(N-1);
that=k*tau/(N-1);

USdot = -bhat*(uI+uA)*uS - that*uQ*uS + delta*uTs;
uEdot =  bhat*(uI+uA)*uS - that*uQ*uE - (alphaA+alphaI)*uE;
uAdot =  alphaA*uE - gammaA*uA - that*uQ*uA;
uIdot =  alphaI*uE - gammaI*uI - (that*uQ+sigma)*uI;
uRdot =  gammaA*uA + gammaI*uI - that*uQ*uR + delta*uTr;
uTsdot=  that*uQ*uS - delta*uTs;
uTrdot=  that*uQ*uR - delta*uTr;
uTqdot=  that*(uE+uA+uI)*uQ + sigma*uI - delta*uTq;
uQdot =  delta*uTq - gammaQ*uQ;
uRqdot=  gammaQ*uQ;

udot=[USdot uEdot uAdot uIdot uRdot uTsdot uTrdot uTqdot uQdot uRqdot]';

