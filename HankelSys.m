classdef HankelSys < handle & matlab.mixin.Copyable
    properties
        sys     %state space model
        n       %system order
        history %Past I/O data
        H       %Hankel matrix
        par     %parameters for a Hankel Matrix
        x0      %initial condition
    end
    methods

        function initialize(obj)
            if isempty(obj.sys)
                n1=obj.n;
            else
                n1=size(obj.sys.A,1);
                obj.n=n1;
            end
            obj.x0=zeros(n1,1);
        end

        function y = stepForward(obj, u,noise)
        y = obj.sys.C * obj.x0 + obj.sys.D * u+noise;
        x = obj.sys.A * obj.x0 + obj.sys.B * u;
        obj.x0=x;
        obj.history = [obj.history, [u;y]];
            if size(obj.history,2) > obj.par.Td
                obj.history(:,1) = [];  
            end
        end

        function makeHankel(obj)
            obj.H=creatHankel(obj.history,obj.par);
        end

        function out=iscomplete(obj)
            r=rank([obj.H.U;obj.H.Y])-(obj.par.m * obj.par.L+obj.n);
            if r==0
                out=["Positive",num2str(rank([obj.H.U;obj.H.Y]))];
            else
                out="Negative";
            end
        end

        function persistent_excitation(obj)
            Td=obj.par.Td;
            m=obj.par.m;
            u=randn(m,Td);
            xini=randn(obj.n,1);
            yout=lsim(obj.sys,u,1:Td,xini);
            obj.history=[u;yout'];
        end

    end

end
