clear;clc;
%设定初值%
mothernum=200;childnum=5;
p=0.6;q=0.4;
%uo=1;ao=0.7;do=0.4;u1=0.1;u2=0.05;v1=0.05;v2=0.10;
%imprint=0.15;
uo=1;ao=0.7;do=0.4;u1=0;u2=0;v1=0;v2=0;
imprint=0;
mothfreq=[p^2 2*p*q q^2];
childfreq=[p q 0 0;0.5*p 0.5*q 0.5*p 0.5*q;0 0 p q];

Uo=[uo+ao+u1 uo+do+imprint+v1 0 0;
    uo+ao-u1 uo+do+imprint-v1 uo+do-imprint-v2 uo-ao-u2;
    0 0 uo+do-imprint+v2 uo-ao+u2];
clear um am dm u0 ao do imprint u1 u2 v1 v2 uo

%产生数据%
t=1;
while t<=100
mothgenotype=randsrc(mothernum,1,[1:3;mothfreq]);
    for i=1:mothernum
        genotype=mothgenotype(i);
        chigenotype=randsrc(1,childnum,[1:4;childfreq(genotype,:)]);
        childgenotype(i,1:childnum)=chigenotype;
    end

for i=1:mothernum
    for j=1:childnum
       simupheno(i,j)=Uo(mothgenotype(i,1),childgenotype(i,j));
    end
end    
for i=1:childnum
fchildgenotype(:,i)=childgenotype(:,i)+10*mothgenotype;
end
simupheno=reshape(simupheno,1,[]);fchildgenotype=reshape(fchildgenotype,1,[]);
%variance=var(simupheno(fchildgenotype==12|fchildgenotype==22|fchildgenotype==23|fchildgenotype==32));
variance=0.04;
clear simupheno fchildgenotype

heritabilitylevel=0.05;
for i=1:mothernum
    for j=1:childnum
      childphenotype(i,j)=normrnd(Uo(mothgenotype(i,1),childgenotype(i,j)),(variance/heritabilitylevel)^0.5);
    end
end

for i=1:childnum
    childgenotype(childgenotype(:,i)==3,i)=2;
    childgenotype(childgenotype(:,i)==4,i)=3;
end

for i=1:mothernum
    data(2*i-1,1,t)=mothgenotype(i,1);
    data(2*i-1,2:1+childnum,t)=childgenotype(i,:);
    data(2*i,2:1+childnum,t)=childphenotype(i,:);
end
t=t+1
end
save('2new2005unifyi0.05.mat', 'data');