function [H]=HankelMatrix(w,L,Tini)
%L: length of the block hankel matrix
%Tini: Partition
%L-Tini=Tf
H.U=hankel(w(1:L),w(L:end));
H.Up=H.U(1:Tini,:);
H.Uf=H.U(Tini+1:end,:);
end