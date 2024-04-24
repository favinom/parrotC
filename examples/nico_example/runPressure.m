clear all
close all

for i=1:41
    
    [int_r,int_i]=computePressure(i);
    
    ir(i)=int_r;
    ii(i)=int_i;
    
end

omega=2*pi*10.^(-5:0.25:5);

loglog(omega,ir,omega,ii);
%loglog(omega,ii)