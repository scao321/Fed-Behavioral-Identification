% ================================================================
%  Federated Behavioral Identification via Grassmann Aggregation
%  Reproduce Fig. 4: Federated SysID and Trajectory Simulation
%  Author: Shiang Cao
%  ================================================================
clear
clc
close all
rng(123)
base_model=c2d(rss(4,2,2),1);
m=2;p=2;Td=200; L=40; Tini=10;
par=struct('m',m,'p',p,'Td',Td,'L',L,'Tini',Tini);
%par.m=2;par.p=2; par.Td=100; par.L=20; par.Tini=5;

sys_base=setup(par,base_model);
P=proj([sys_base.H.Up;sys_base.H.Uf;sys_base.H.Yp;sys_base.H.Yf]);

num=100;
for i=1:num
    model1=peturb(base_model, 0.1,0.1);
    sys1=setup(par,model1);
    H=[sys1.H.Up;sys1.H.Uf;sys1.H.Yp;sys1.H.Yf];
    Pe(:,:,i)=proj(H);
    %L_gap(i)=norm(Pe(:,:,i)-P,2);
end

%mean(L_gap)
%%

Pavg=proj(mean(Pe,3),m*L+4);
%norm(P-Pavg,2)

%%
% simulation test
yini=zeros(p*Tini,1);
uini=zeros(m*Tini,1);
rng(123)
%uf=randn(L-Tini,2);
uf=randn(60,2);
Lf=size(uf,1);
%yd=ddsimulate([sys_base.H.Up;sys_base.H.Uf;sys_base.H.Yp;sys_base.H.Yf],par,uf);

yd=lsim(base_model,uf,1:Lf);

yd2=ddsimulate(Pavg,par,uf);

yd3=ddsimulate(Pe(:,:,1),par,uf);
yd4=ddsimulate(Pe(:,:,10),par,uf);



% colororder(tab10)
% subplot(211)
% plot(1:Lf,yd(:,1),"LineWidth",1)
% hold on
% plot(1:Lf,yd2(:,1),LineWidth=1,LineStyle="-")
% plot(1:Lf,yd3(:,1),"LineWidth",1,"LineStyle","--")
% grid on
% grid minor
% subplot(212)
% plot(1:Lf,yd(:,2),"LineWidth",1)
% hold on
% plot(1:Lf,yd2(:,2),LineStyle="-",LineWidth=1)
% plot(1:Lf,yd3(:,2),"LineStyle","--",LineWidth=1)
% plot(1:Lf,yd4(:,2),"LineStyle","-.",LineWidth=1)
% legend('Ture','FedAvg','Similar System1','Similar system 2')
% grid on
% grid minor
% 
% set(gca,'fontsize',15)


%%
for i=1:num
    yout=ddsimulate(Pe(:,:,i),par,uf);
    y1(i,:)=yout(:,1)';
    y2(i,:)=yout(:,2);
end
[N, T] = size(y1);

upper_env = max(y1, [], 1);   
lower_env = min(y1, [], 1);  

t = 1:T;

figure

subplot(211)
hold on
colororder(tab10)
shade=fill([t, fliplr(t)], [upper_env, fliplr(lower_env)], ...
     [0.9,0.9 .9], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
h1=plot(t,yd(:,1),LineWidth=1.5);
h2=plot(t,yd2(:,1),LineWidth=1.5,LineStyle="--");
grid on
grid minor
%xlabel('Time');
ylabel('y_1');
set(gca,'fontsize',15)
legend([shade,h1,h2],{'Envelope of the similar systems','Target system','Federated estimated system'},fontsize=12)


upper_env2 = max(y2, [], 1);   
lower_env2 = min(y2, [], 1); 
subplot(212)
hold on
colororder(tab10)
fill([t, fliplr(t)], [upper_env2, fliplr(lower_env2)], ...
     [0.9,0.9 .9], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
plot(t,yd(:,2),LineWidth=1.5)
plot(t,yd2(:,2),LineWidth=1.5,LineStyle="--")
xlabel('Time');
ylabel('y_2');
set(gca,'fontsize',15)
grid on
grid minor

%%
function y=ddsimulate(P,par,uf) 
% uf [Tf,m]
L_uf=size(uf,1);
m=par.m; p=par.p; L=par.L; Tini=par.Tini;
yini=zeros(p*Tini,1);
uini=zeros(m*Tini,1);
y=[];
for i=1:30:L_uf
    uff=reshape(uf(i:i+29,:),[],1);
g=pinv(P(1:m*L+p*Tini,:))*[uini;uff;yini];
yd=P(m*L+p*Tini+1:end,:)*g;
yini=yd(end-p*Tini+1:end,:);
uini=uff(end-m*Tini+1:end,:);
y=[y;reshape(yd,[],p)];
end

end


function P=proj(H,k)
[U,~,~] = svd(H, 'econ'); 

if nargin<2
    k=rank(H);
end

Uproj = U(:, 1:k);         
P = Uproj * Uproj';   
end

function sys=setup(par,ss_model)
sys=HankelSys;
sys.par=par;
sys.sys=ss_model;
sys.initialize();
sys.persistent_excitation();
sys.makeHankel
end


function ss_model=peturb(base_model, var1,var2)
A=base_model.A; B=base_model.B; C=base_model.C; D=base_model.D;
Ts=base_model.Ts;
dA= var1*randn(size(A)); dB= var2*randn(size(B));
ss_model=ss(A+dA,B+dB,C,D,Ts);
end


