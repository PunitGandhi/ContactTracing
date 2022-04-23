close all;
clear all;


beta=0.05; %probability connection w/infected leads to exposure

alphaI=0.6; %probability symptomatic
alphaA=0.1; %probability asymptomatic

gammaA=0.1; %probability of recovery Asymptomatic
gammaI=0.1; %probability of recovery Infected
gammaQ=0.5; %probability of quarantine (and tracing) complete

sigma=0.3; %probability infected will decide to get tested

tau=gammaQ; %probability connection of tested, fix to gammaQ
delta=0.3; %probability test is complete

dt=0.1; %time step in days


k0=20; %use to set beta and tau
beta0=0.05;
tau0=gammaQ;

NN=1000; %number of nodes
WSp=0.1; %0.01 0.1

Nkk=1;
kvals=10;

NT=round(200/dt); %number of time steps not including initial value

Nll=200;
MaxParForN=100;

OutbreakThreshold=20;
k=kvals;
fnm=['networkODEcompareWSp' num2str(100*WSp) 'Deg' num2str(2*k)...
        'tau' num2str(tau0*100) 'sigma' num2str(sigma*100)...
        'Nodes' num2str(NN) 'runs' num2str(Nll)];



ItotQ=zeros(Nkk,Nll);
ItotR=zeros(Nkk,Nll);
Iend=zeros(Nkk,Nll);
Ipeak=zeros(Nkk,Nll);
Tpeak=zeros(Nkk,Nll);
Npeak=zeros(Nkk,Nll);
Tlength=zeros(Nkk,Nll);
totTests=zeros(Nkk,Nll);
selfTests=zeros(Nkk,Nll);
tic

%for (ll=1:Nll) %number of different networks to try
    
for ll=1:Nll %number of different networks to try
    ItotQnet=zeros(1,Nkk);
    ItotRnet=zeros(1,Nkk);
    Iendnet=zeros(1,Nkk);
    Ipeaknet=zeros(1,Nkk);
    Tpeaknet=zeros(1,Nkk);
    Npeaknet=zeros(1,Nkk);
    Tlengthnet=zeros(1,Nkk);
    totTestsnet=zeros(1,Nkk);
    selfTestsnet=zeros(1,Nkk);
    
    for kk=1:length(kvals)
        
        k=kvals(kk); %2k is avg degree
        
        
        beta=beta0*k0/k;
        tau=tau0*k0/k;
        %G = ErdosRenyi(NN,k); %Random network
        G = WattsStrogatz(NN,k,WSp); %small world network
        %G = RandRegGraph(NN,k); %Random network
        %G=graph(ones(NN,NN)-diag(ones(NN,1)));
        %figure;
        %plot(G,'NodeColor','k','EdgeAlpha',0.1);
        Madj=adjacency(G);
        
        
        
        
        
        
        %%initial conditiontspan
        SSt=true(NN,1);
        EEt=false(NN,1);
        IIt=false(NN,1);
        AAt=false(NN,1);
        RRt=false(NN,1);
        QQt=false(NN,1);
        RQt=false(NN,1);
        
        TSt=false(NN,1);
        TPt=false(NN,1);
        TRt=false(NN,1);
        
        tt=0;
        tlast=0;
        
        %infect random node(s)
        ninfect=1; % could be less because of repeats
        n0infect=unique(ceil(rand(ninfect,1)*(NN)));
        EEt(n0infect,1)=true;
        SSt(n0infect,1)=false;
        
        %save init cond
        SS=SSt;
        EE=EEt;
        II=IIt;
        AA=AAt;
        RR=RRt;
        QQ=QQt;
        RQ=RQt;
        
        TS=TSt;
        TP=TPt;
        TR=TRt;
        
        ntcount=-1;
        testcount=0;
        selftestcount=0;
        Iflag=1;
        while Iflag>0
            
            ntcount=ntcount+1;
            %States 1 if node in this state 0 otherwise.
            SStmp=false(NN,NT);
            EEtmp=false(NN,NT); %currently not infectious
            IItmp=false(NN,NT);
            AAtmp=false(NN,NT);
            RRtmp=false(NN,NT);
            QQtmp=false(NN,NT);
            RQtmp=false(NN,NT);
            
            TStmp=false(NN,NT);
            TRtmp=false(NN,NT);
            TPtmp=false(NN,NT);
            
            tttmp=zeros(1,NT);
            
            
            %for nt=(2:(NT+1)) + (NT)*ntcount
            for nt=(1:(NT))
                rn=rand(NN,1); %use to decide transitions
                
                %probability of susceptible not getting exposed
                %Ps2s=((1-beta).^(Madj*(IIt+AAt)));
                %Ps2s=1-beta*(Madj*(IIt+AAt));
                %indSE=(rn>Ps2s) & SSt;
                %probability of susceptible getting exposed
                Ps2e=beta*(Madj*(IIt+AAt))*dt;
                indSE=(rn<Ps2e) & SSt;
                
                %probability of exposed developing symptoms;
                Pe2i=alphaI*EEt*dt;
                indEI=rn<Pe2i; %& EEt;
                %probability of exposed becoming asymptomatic
                Pe2a=Pe2i+alphaA*EEt*dt;
                indEA=(rn<Pe2a) & (~indEI); %& EEt;
                
                %probability of infected recovered
                Pi2r=gammaI*IIt*dt;
                indIR=(rn<Pi2r); %& IIt;
                %probability of infected self-quarentine
                %Pi2q=Pi2r+sigma*IIt;
                %indIQ=rn<Pi2q & (rn>Pi2r) & IIt;
                
                %probability of asymptomatic recovered
                Pa2r=gammaA*AAt*dt;
                indAR=(rn<Pa2r); % & AAt;
                
                %probability of quarantined complete
                Pq2r=gammaQ*QQt*dt;
                indQR=(rn<Pq2r); % & QQt;
                
                %probability of not getting tested from tracing
                %Px2t=((1-tau).^(Madj*(QQt)));
                %probability of getting tested from tracing
                Px2t=tau*(Madj*(QQt))*dt;
                Pt2x=delta*dt;
                
                %probability of symptomatic choosing to get tested
                Pi2t=sigma*IIt*dt;
                
                indST=( (rn)<(Ps2e + Px2t) )  & (~indSE) & SSt;
                indTS=(rn< Pt2x) & TSt;
                indET=(rn<(Pe2a+Px2t)) & (~(indEI|indEA)) & EEt;
                indIT=(rn<(Pi2r+Px2t+Pi2t)) & (~indIR) & IIt;
                indAT=(rn<(Pa2r+Px2t)) & (~indAR) & AAt;
                indTQ=(rn<Pt2x) & TPt;
                indRT=(rn<Px2t) & RRt;
                indTR=(rn<Pt2x) & TRt;
                
                %count up tests used
                testcount=testcount+sum(indST+indET+indIT+indAT+indRT);
                %sum(QQt)
                stst=indIT & (rn>(Pi2r+Px2t));
                selftestcount=selftestcount+sum(stst);
                %Update S-> E
                EEt(indSE)=true;
                SSt(indSE)=false;
                
                %update E-> I
                IIt(indEI)=true;
                EEt(indEI)=false;
                %update E->A
                AAt(indEA)=true;
                EEt(indEA)=false;
                
                %update I-> R
                RRt(indIR)=true;
                IIt(indIR)=false;
                %update I-> Q
                %QQt(indIQ)=true;
                %IIt(indIQ)=false;
                
                %update A->R
                RRt(indAR)=true;
                AAt(indAR)=false;
                
                %update Q->R
                RQt(indQR)=true;
                QQt(indQR)=false;
                
                %Update S->T
                TSt(indST)=true;
                SSt(indST)=false;
                
                %Update T->S
                TSt(indTS)=false;
                SSt(indTS)=true;
                
                %Update E->T
                TPt(indET)=true;
                EEt(indET)=false;
                %Update I->T
                TPt(indIT)=true;
                IIt(indIT)=false;
                %Update A->T
                TPt(indAT)=true;
                AAt(indAT)=false;
                
                %Update T->Q
                TPt(indTQ)=false;
                QQt(indTQ)=true;
                
                %Update R->T
                TRt(indRT)=true;
                RRt(indRT)=false;
                
                %Update T->R
                TRt(indTR)=false;
                RRt(indTR)=true;
                
                
                
                
                
                
                
                %store outputs
                SStmp(:,nt)=SSt;
                EEtmp(:,nt)=EEt;
                AAtmp(:,nt)=AAt;
                IItmp(:,nt)=IIt;
                RRtmp(:,nt)=RRt;
                QQtmp(:,nt)=QQt;
                RQtmp(:,nt)=RQt;
                
                TStmp(:,nt)=TSt;
                TRtmp(:,nt)=TRt;
                TPtmp(:,nt)=TPt;
                
                tlast=tlast+dt;
                tttmp(nt)=tlast;
                
                
            end
            
            %if Iflag>0
            SS=[SS, SStmp];
            EE=[EE, EEtmp];
            II=[II, IItmp];
            AA=[AA, AAtmp];
            RR=[RR, RRtmp];
            QQ=[QQ, QQtmp];
            RQ=[RQ, RQtmp];
            
            TS=[TS,TStmp];
            TR=[TR,TRtmp];
            TP=[TP,TPtmp];
            
            tt=[tt,tttmp];
            %end
            
            
            Iflag=sum(TPt+IIt+AAt+QQt+EEt);
            Iendnet(kk)=Iflag; %total number with active infection at end of simulation
            %%%%%%%%
            Iflag=0; %stop simulation after NT steps
            %%%%%%%%
        end
        ntt=length(tt);
                %Store time series data
        nukS(1:ntt,kk)=sum(SS,1);
        nukE(1:ntt,kk)=sum(EE,1);
        nukA(1:ntt,kk)=sum(AA,1);
        nukI(1:ntt,kk)=sum(II,1);
        nukR(1:ntt,kk)=sum(RR,1);
        nukQ(1:ntt,kk)=sum(QQ,1);
        nukRq(1:ntt,kk)=sum(RQ,1);
        nukTs(1:ntt,kk)=sum(TS,1);
        nukTr(1:ntt,kk)=sum(TR,1);
        nukTp(1:ntt,kk)=sum(TP,1);
        Tkend(kk)=tt(end);
        
        ItotQnet(kk)=sum(RQt); %total number infected is total number recovered
        ItotRnet(kk)=sum(RRt); %total number infected is total number recovered
        
        totI=sum(TP+II+AA+QQ+EE);
        [maxtmp, indtmp]= max(totI);
        [pks,locs,w]=findpeaks(totI,'MinPeakProminence',OutbreakThreshold);
        
        Npeaknet(kk)=length(pks);
        Ipeaknet(kk)=maxtmp; %peak number of people who infected at a given time
        
        Tpeaknet(kk)=tt(indtmp);
        
        Tlengthnet(kk)=tt(find(totI==0,1,'first'));
        totTestsnet(kk)=testcount;
        selfTestsnet(kk)=selftestcount;
        
    end
    

            %Store time series data
        nuS(:,:,ll)=nukS;
        nuE(:,:,ll)=nukE;
        nuA(:,:,ll)=nukA;
        nuI(:,:,ll)=nukI;
        nuR(:,:,ll)=nukR;
        nuQ(:,:,ll)=nukQ;
        nuRq(:,:,ll)=nukRq;
        nuTs(:,:,ll)=nukTs;
        nuTr(:,:,ll)=nukTr;
        nuTp(:,:,ll)=nukTp;
        Tend(:,:,ll)=Tkend;
        
        
        
    
    Iend(:,ll)=Iendnet; %total number with active infection at end of simulation
    
    
    
    
    ItotQ(:,ll)=ItotQnet; %total number infected is total number recovered
    ItotR(:,ll)=ItotRnet; %total number infected is total number recovered
    
    Npeak(:,ll)=Npeaknet;
    Ipeak(:,ll)=Ipeaknet; %peak number of people who infected at a given time
    
    Tpeak(:,ll)=Tpeaknet;
    
    Tlength(:,ll)=Tlengthnet;
    totTests(:,ll)=totTestsnet;
    selfTests(:,ll)=selfTestsnet;
end
toc


%%
%% Parameters
par.k=2*k;
par.beta=beta;%% Plot results
par.N=NN;
par.tau=tau;
par.delta=delta;
par.alphaA=alphaA;
par.alphaI=alphaI;
par.gammaA=gammaA;
par.gammaI=gammaI;
par.sigma=sigma;
par.gammaQ=gammaQ;


%% initial Condition
uE0=ninfect;

uS0=par.N-uE0;
%uE0=;
uA0=0;
uI0=0;
uR0=0;
uTs0=0;
uTr0=0;
uTq0=0;
uQ0=0;
uRq0=0;

u0=[uS0 uE0 uA0 uI0 uR0 uTs0 uTr0 uTq0 uQ0 uRq0];



%% Run ODE Simulation
tspan=(0:NT)*dt;

[t,u] = ode15s(@ContactTracingODEmodel, tspan, u0, [], par);



%% Plot results

uS=u(:,1);
uE=u(:,2);
uA=u(:,3);
uI=u(:,4);
uR=u(:,5);
uTs=u(:,6);
uTr=u(:,7);
uS=u(:,1);
uE=u(:,2);
uA=u(:,3);
uI=u(:,4);
uR=u(:,5);
uTs=u(:,6);
uTr=u(:,7);
uTq=u(:,8);
uQ=u(:,9);
uRq=u(:,10);



uTq=u(:,8);
uQ=u(:,9);
uRq=u(:,10);



save([fnm 'simdata.mat']);


pltg=[0.6 0.6 0.6 0.1];
pltYmax=[1000, 100, 100, 750, 100];
indlvals=1:Nll;
indkvals=1:Nkk;
fntsz=16;

%%%CREAT FIGURE
figure('Position',[100 100, 400, 800]);
tspan=(0:NT)*dt;
indt=1:length(tspan);


uS=u(:,1);
uE=u(:,2);
uA=u(:,3);
uI=u(:,4);
uR=u(:,5);
uTs=u(:,6);
uTr=u(:,7);
uTq=u(:,8);
uQ=u(:,9);
uRq=u(:,10);






%%S
subplot(5,1,1)
hold on;box on;
set(gca,'FontSize',fntsz,'defaultTextInterpreter','latex')

for kk=indkvals
    for ll=indlvals
        plot(tspan,nuS(indt,kk,ll),'Color', pltg);
    end
    %plot(tspan,median(nuS(indt,kk,:),3),'k--')
end
%tspan=1:NT;
plot(t,uS,'LineWidth',2,'Color','blue');
%plot(t,duS,'LineWidth',2,'Color','red','LineStyle',':');
ylabel('$S$')
ylim([0,pltYmax(1)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,1,2)
hold on;box on;
set(gca,'FontSize',fntsz,'defaultTextInterpreter','latex')
for kk=indkvals
    for ll=indlvals
        plot(tspan,nuE(indt,kk,ll),'Color', pltg);
    end
    %plot(tspan,median(nuE(indt,kk,:),3),'k--')
end
%tspan=1:NT;
plot(t,uE,'LineWidth',2,'Color','blue');
%plot(t,duE,'LineWidth',2,'Color','red','LineStyle',':');
ylabel('$E$')
ylim([0,pltYmax(2)])
%%%%%%%%%%box on;
set(gca,'FontSize',fntsz,'defaultTextInterpreter','latex')%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,1,3)
hold on;box on;
set(gca,'FontSize',fntsz,'defaultTextInterpreter','latex')
for kk=indkvals
    for ll=indlvals
        plot(tspan,nuA(indt,kk,ll)+nuI(indt,kk,ll),'Color', pltg);
    end
    %plot(tspan,median(nuI(indt,:,kk)+nuA(indt,kk,:),3),'k--')
end
%tspan=1:NT;
plot(t,uA+uI,'LineWidth',2,'Color','blue');
%plot(t,duA+duI,'LineWidth',2,'Color','red','LineStyle',':');
ylabel('$I+A$')
ylim([0,pltYmax(3)])
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,1,4)
hold on;box on;
set(gca,'FontSize',fntsz,'defaultTextInterpreter','latex')
for kk=indkvals
    for ll=indlvals
        plot(tspan,nuTs(indt,kk,ll)+nuTr(indt,kk,ll)+nuTp(indt,kk,ll),'Color', pltg);
    end
    %plot(tspan,median(nuTs(indt,kk,:)+nuTr(indt,kk,:)+nuTp(indt,kk,:),3),'k--')
end
%tspan=1:NT;
plot(t,uTs+uTr+uTq,'LineWidth',2,'Color','blue');
%plot(t,duTs+duTr+duTq,'LineWidth',2,'Color','red','LineStyle',':');
ylabel('$T$')
ylim([0,pltYmax(4)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,1,5)
hold on;box on;
set(gca,'FontSize',fntsz,'defaultTextInterpreter','latex')
for kk=indkvals
    for ll=indlvals
        plot(tspan,nuQ(indt,kk,ll),'Color', pltg);
    end
    %plot(tspan,median(nuQ(indt,kk,:),3),'k--')
end
%tspan=1:NT;
plot(t,uQ,'LineWidth',2,'Color','blue');
%plot(t,duQ,'LineWidth',2,'Color','red','LineStyle',':');
ylabel('$Q$')
ylim([0,pltYmax(5)])
xlabel('$t$ (days)')
print([fnm 'tseries.png'],'-dpng')


