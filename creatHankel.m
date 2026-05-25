function [sys]=creatHankel(w,par)
% w trajectory [U;Y],each row is the trajectory
% par parameters used to specify the Hankel Matrix
if isfield(par, 'L')
    L=par.L;
else
    L=par.Tini+par.Tf;
end

if ~isfield(par,'Tini')
    Tini=[];
else
    Tini=par.Tini;
end



m=par.m; %number of inputs
[nor,~]=size(w);


U=[];
Up=[];
Uf=[];
for i=1:m
    H=HankelMatrix(w(i,:),L,Tini);
    sys.(sprintf('U%d', i))=H;  
    U=[U;H.U];
    Up=[Up;H.Up]; %[u1_p;u2_p;...]
    Uf=[Uf;H.Uf]; %[u1_f;u2_fl...
end

Y=[];
Yp=[];
Yf=[];
for i = 1:nor-m
    H=HankelMatrix(w(i+m,:),L,Tini);
    sys.(sprintf('Y%d', i)) =H;
    Y=[Y;H.U];
    Yp=[Yp;H.Up];
    Yf=[Yf;H.Uf];
end

sys.U=U;
sys.Y=Y;
sys.Up=Up;
sys.Uf=Uf;
sys.Yf=Yf;
sys.Yp=Yp;

end