% ================================================================
%  Federated Behavioral Identification via Grassmann Aggregation
%  Reproduce Fig. 3: Structured perturbation with varying number of local systems
%  Author: Shiang Cao
%  ================================================================
clear
clc
rng(123)
base_model=c2d(rss(4,2,2),1);
m=2;p=2;Td=200; L=40; Tini=10;
par=struct('m',m,'p',p,'Td',Td,'L',L,'Tini',Tini);
%par.m=2;par.p=2; par.Td=100; par.L=20; par.Tini=5;

sys_base=setup(par,base_model);
P=proj([sys_base.H.Up;sys_base.H.Uf;sys_base.H.Yp;sys_base.H.Yf]);

totalnum=[5:5:100];
for j=1:length(totalnum)
num=totalnum(j);
Pe=[];
for i=1:num
    model1=peturb(base_model, 0.1,0.1);
    sys1=setup(par,model1);
    H=[sys1.H.Up;sys1.H.Uf;sys1.H.Yp;sys1.H.Yf];
    Pe(:,:,i)=proj(H);
    L_gap(i)=norm(Pe(:,:,i)-P,2);
end

d1(j)=mean(L_gap);

Pavg=proj(mean(Pe,3),m*L+4);
d2(j)=norm(P-Pavg,2);

j
end
%%
plot(totalnum,d1,Marker="*",LineWidth=1.5)
hold on
plot(totalnum,d2,Marker="o",LineWidth=1.5)
grid on
grid minor
xlabel('Number of local systems')
ylabel('L-gap')
set(gca,'fontsize',15)
legend('Mean of the local systems','Federated estimate')


%%
function y=ddsimulate(P,par,uf) 
% uf [Tf,m]
m=par.m; p=par.p; L=par.L; Tini=par.Tini;
yini=zeros(p*Tini,1);
uini=zeros(m*Tini,1);
g=pinv(P(1:m*L+p*Tini,:))*[uini;reshape(uf,[],1);yini];
yd=P(m*L+p*Tini+1:end,:)*g;
y=reshape(yd,[],p);
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


