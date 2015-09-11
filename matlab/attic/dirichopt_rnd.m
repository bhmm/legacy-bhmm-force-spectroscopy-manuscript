function ddraw = dirichopt_rnd(a)
% PURPOSE: a matrix of random draws from the Dirichlet distribution
%---------------------------------------------------

kdim=size(a,1);
a1=zeros(kdim,1);
for i = 1:kdim
%    a1(i,1)=gamm_rnd(1,1,a(i,1),1);
    %a1(i,1)=gamm_rnd(a(i,1));
  m=a(i);
  if m<=1
    % Use RGS algorithm by Best, p. 426
    c=1/m; 
    t=0.07+0.75*sqrt(1-m);
    b=1+exp(-t)*m/t;
    %for i1=1:nrow
    %  for i2=1:ncol
         accept=0;
         while accept==0
            u=rand; w=rand; v=b*u;
            if v<=1
               x=t*(v^c);
               accept=((w<=((2-x)/(2+x))) | (w<=exp(-x)));
            else
               x=-log(c*t*(b-v));
               y=x/t;
               accept=(((w*(m+y-m*y))<=1) | (w<=(y^(m-1))));
            end
         end
    %     gb(i1,i2)=x;
         a1(i)=x;
    %  end
    %end
  else
    % Use Best's rejection algorithm XG, p. 410
    b=m-1;
    c=3*m-0.75;
    %for i1=1:nrow
    %  for i2=1:ncol
         accept=0;
         while accept==0
            u=rand;  v=rand;
            w=u*(1-u);  y=sqrt(c/w)*(u-0.5);
            x=b+y;
            if x >= 0
               z=64*(w^3)*v*v;
               accept=(z<=(1-2*y*y/x)) ...
                      | (log(z)<=(2*(b*log(x/b)-y)));
            end
         end
    %     gb(i1,i2)=x;
         a1(i)=x;
    %  end
    %end
  end
  %gb=gb/k;    
  %a1(i,1)=gb;
end
ddraw=a1./sum(a1);
