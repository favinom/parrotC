close all
clear all

a=ls;

k = strfind(a,'cleaned');

[m,n]=size(k);
m=m*n;

if (m<1)

    disp('here')
    
!./cleanFile.sh level2.csv
!./cleanFile.sh level3.csv
!./cleanFile.sh level4.csv
!./cleanFile.sh level5.csv
!./cleanFile.sh level6.csv
!./cleanFile.sh level7.csv

%!./cleanFile.sh direct.csv

!touch cleaned

end

R=csvread('direct.csv');
A{1}=csvread('level2.csv');
A{2}=csvread('level3.csv');
A{3}=csvread('level4.csv');
A{4}=csvread('level5.csv');
A{5}=csvread('level6.csv');
A{6}=csvread('level7.csv');

plot(A{1}(:,1),R(:,6))
hold on
plot(A{1}(:,1),A{1}(:,6))
plot(A{1}(:,1),A{2}(:,6))
plot(A{1}(:,1),A{3}(:,6))
plot(A{1}(:,1),A{4}(:,6))
plot(A{1}(:,1),A{5}(:,6))
plot(A{1}(:,1),A{6}(:,6))


legend('direct','level2','level3','level4','level5','level6','level7')