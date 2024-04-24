clear all
close all

a=dir('newoutput_*b.csv');

datab=[];
refb=[];

for i=1:length(a)
    
    str=a(i).name;
    newStr = split(str,'_');
    
    refb=[refb; str2num(newStr{2}) str2num(newStr{3})];
    datab{i}=csvread(str);
    
end

a=dir('newoutput_*f.csv');

dataf=[];
reff=[];

for i=1:length(a)
    
    str=a(i).name;
    newStr = split(str,'_');
    
    reff=[refb; str2num(newStr{2}) str2num(newStr{3})];
    dataf{i}=csvread(str);
    
end

freq=dataf{1}(:,1);

subplot(1,2,1)
plot(freq,dataf{1}(:,2))
hold on
plot(freq,dataf{2}(:,2))
plot(freq,dataf{3}(:,2))
plot(freq,dataf{4}(:,2))
plot(freq,dataf{5}(:,2))

plot(freq,dataf{6}(:,2),'-')

subplot(1,2,2)
plot(freq,datab{1}(:,2),'b')
hold on
plot(freq,datab{2}(:,2),'b')
plot(freq,datab{3}(:,2),'b')
plot(freq,datab{4}(:,2),'b')
plot(freq,datab{5}(:,2),'b')

plot(freq,datab{6}(:,2),'r-')
plot(freq,datab{7}(:,2),'r-')
plot(freq,datab{8}(:,2),'r-')
plot(freq,datab{9}(:,2),'r-')
plot(freq,datab{10}(:,2),'r-')

figure(2)

plot(freq,dataf{15}(:,2),freq,datab{15}(:,2))

return

return
