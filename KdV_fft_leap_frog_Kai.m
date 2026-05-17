%% KdV_fft_leap_frog_Kai.m
% solve the Kdv equation
% use the mass conservative scheme
% fft in space 1/p*(u^p)_x
% 2nd order leap-frog scheme in time

clear;
% clc
tic;

% L=64*pi*1;
% N=2^13;

L=10*pi;
N=2^10;

p=2;
sigma=(p-1)/2;
x=linspace(-L,L,N+1)';
x=x(1:N);
dx=x(2)-x(1);
% for jjjj=0:0
% dt0=2e-2/2^jjjj;
dt0=1e-2/2^0;
% dt0=1e-2/2^8;
t=0;
tmax=50;

tplot=t;
Ntplot=round(1/dt0); % save the time at every 100 iterations
% Ntplot=round(10);
Nt=1e5; % simulate 1e5 time steps
k=[0:N/2-1 0 -N/2+1:-1]';
k1=1i*pi/L*k; % 1st derivative on frequency space
k2=k1.^2;
k3=k1.*k2;
k4=k2.^2;
set(0,'defaultAxesFontSize',16)
%% set up the initial conditions and rescaling
Qf=@(r)(p*(p+1)/2)^(1/(p-1))*(sech((p-1)/2*r)).^(2/(p-1));
Qcf=@(r,c) c^(1/(p-1))*Qf(c^.5*r); 

xs1=-0;
xs2=+15;
c1=1;
c2=0.5;
QQ1=Qcf(x-xs1,c1)*1;
QQ2=Qcf(x-xs2,c2);
u0=1*QQ1-0*QQ2;

vQx=mod(x-c1*t+L+xs1,2*L)-L;
gamma=c1;
uexactf=@(t,vQx) c1^(1/(p-1))*Qf(c1^.5*vQx); 
uexact=u0;
%% 
uh0=fft(u0);
u_x=real(ifft(uh0.*k1));
%% find the initial mass and energy
energy=real(0.5*sum(-k2.*uh0.*conj(uh0))*2*L/N^2-1/(p*(p+1))*sum(u0.^(p+1))*2*L/N);
mass=sum(abs(u0).^2)*2*L/N;

EM=0; % error of mass
EE=0; % error of energy
ue=0; % error of the solution

Q=Qcf(x,1);
Qh=fft(Q);
MQ=sum(Q.^2)*2*L/N;
EQ=.5*sum(-k2.*Qh.*conj(Qh))*2*L/N^2-1/(p*(p+1))*sum(Q.^(p+1))*2*L/N;
Eu=energy;
%% save data
u00=u0;
uh00=uh0;
U=u00;
UE=u00; % the exact solution
% A00=norm(u00,inf); % initial height
[A00,Nc]=max(abs(u00));
dt=dt0; % the time step
AA=A00;
xxc=x(Nc);
%% u1 obtained by IRK2
iterf=0;
errorf=1;
tolf=1e-12;
u05=u0;
while errorf>tolf;
    Nu05=u05.^p;
    Nu05h=fft(Nu05);
    RHSh=uh0-dt/2*1/p*k1.*Nu05h;
    uh051=RHSh./(1+dt/2*k3);
    u051=real(ifft(uh051));

    errorf=max(abs(u051-u05));

    u05=u051;
    uh05=uh051;
    iterf=iterf+1;
    if iterf>=100;
        disp('fixed point does not converge')
        break
    end
end
u1=2*u05-u0;
uh1=2*uh05-uh0;
%% start the time iteration
% vQx1=mod(x-c1*dt+L+xs1,2*L)-L;
% u1=uexactf(dt,vQx1);
% uh1=fft(u1);

iteration=0;
t=0;
while t<tmax;
    Nu=u1.^p;
    Nuh=fft(Nu);

    RHSh=1/(2*dt)*uh0-1/2*k3.*uh0-1/p*k1.*Nuh; % 2nd order
    uh2=RHSh./(1/(2*dt)+1/2*k3); %
    u2=real(ifft(uh2));
%% update the information for the next use
    u0=u1; uh0=uh1;
    u1=u2; uh1=uh2;

    t=t+dt;
    iteration=iteration+1;
    A=norm(u0,inf);

    if mod(iteration,Ntplot)==0
        u_x=real(ifft(uh0.*k1));
        energy(end+1)=real(0.5*sum(-k2.*uh0.*conj(uh0))*2*L/N^2-1/(p*(p+1))*sum(u0.^(p+1))*2*L/N);
        mass(end+1)=sum(abs(uh0).^2)*2*L/N^2;


        EM(end+1)=abs(mass(end)-mass(1));
        EE(end+1)=abs(energy(end)-energy(1));
        
        Eu(end+1)=sum((.5*abs(u_x).^2-1/(p*(p+1))*(u0).^(p+1)))*2*L/N;
        tplot(end+1)=t; % time coordinate
        U(:,end+1)=u0;
        
        vQx=mod(x-c1*t+L+xs1,2*L)-L;
        uexact=uexactf(t,vQx);
        
        UE(:,end+1)=uexact;
        ue(end+1)=max(abs(u0-UE(:,end))); % error at time t
        
        AA(end+1)=A;
        xxc(end+1)=x(Nc);
        
        [umax,Nc]=max(u0);
%          plot(x,u0,'-',x,u00,'--','linewidth',2);grid on
         plot(x,u0,x,uexact,'--','linewidth',2);grid on
%          xlim([-25+x(Nc),25+x(Nc)]); 
%          xlim([-50,50]); 
         legend('u(x,t)','exact')
         xlim([-inf,inf])
         ylim([-0,3])
         title(['t=',num2str(t)])
         drawnow;
    end
    
    
    if A>A00*5;
        disp('solution blows up in finite time')
        break
    end    

end
%% display the error
disp(['errorM=', num2str(EM(end))]);
disp(['errorE=', num2str(EE(end))]);
disp(['dt=',num2str(dt0),'; error=',num2str(max(ue))])
disp('**********************************************')
%% plot
figure
plot(tplot,ue,'linewidth',2); grid on
xlim([-inf,inf]); xlabel('t');
legend('L^{\infty} error')

figure % (2)
plot(tplot,EM,tplot,EE,'--','linewidth',2); grid on
legend('error of mass','error of energy')
% set(gca, 'YScale', 'log')

%% save data
totaltime=toc;

toc;



