function [int_r,int_i]=computePressure(input)

% this is kappa over eta
kb=1e-16;
kf=1e-7;

filename=['pressure_',num2str(input),'.csv'];

A=csvread(filename,1);

x=A(:,1);
y=A(:,2);
pi=A(:,4);
pr=A(:,5);

indices=find( abs(x) < 1e-15);

y0=y(indices);

[y0,i]=sort(y0);
indices=indices(i);

counter=1;

for i=1:length(indices)
    
    ii=indices(i);
    
    indicesy=find( abs(y-y(ii)) < 1e-14);
    
    xy=x(indicesy);
    
    if ( 50-0.005/2<=y(ii) && y(ii)<=50+0.005/2 )
        k=kf;
    else
        k=kb;
    end
    
    [xy,i]=sort(xy);
    indicesy=indicesy(i);
    
    ur(counter)=-k*(pr(indicesy(2))-pr(indicesy(1)))/(x(indicesy(2)) - x(indicesy(1)) );
    ui(counter)=-k*(pi(indicesy(2))-pi(indicesy(1)))/(x(indicesy(2)) - x(indicesy(1)) );
    counter=counter+1;
    
end

ur=ur';
ui=ui';

h=diff(y0);

oh=[0;h];
ho=[h;0];

M=spdiags([1/6*ho 1/3*(ho+oh) 1/6*oh], [-1 0 1],length(ur),length(ur));

int_r=sum(M*ur);
int_i=sum(M*ui);
