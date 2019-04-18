clc
clear
cd ('F:\MATLAB\work\');
load('4001unify0.2.mat','data');
mothernum=400;childnum=1;t=1;
for t=1:5
handledata=data(:,:,t);     
for i=1:mothernum
mothgenotype(i,1)=handledata(2*i-1,1);
childgenotype(i,1:childnum)=handledata(2*i-1,2:childnum+1);
childphenotype(i,1:childnum)=handledata(2*i,2:childnum+1);
end
for i=1:childnum
childgenotype(:,i)=10*mothgenotype+childgenotype(:,i);
end
chidata(:,1)=reshape(childgenotype',[],1);
orimn=tabulate(mothgenotype(:));
mn=zeros(3,1);mn(1:numel(orimn(:,1)),1)=orimn(:,2);
oricn=tabulate(chidata(:,1));
cn=zeros(3,3);cn(1,1)=oricn(11,2);cn(1,2)=oricn(12,2);
cn(2,1)=oricn(21,2);cn(2,2)=oricn(22,2);cn(2,3)=oricn(23,2);
cn(3,2)=oricn(32,2);cn(3,3)=oricn(33,2);
chidata(:,2)=reshape(childphenotype',[],1);

%estimate p q%
ep(t,1)=(mn(1)*2+mn(2)+cn(1,1)+cn(2,1)+cn(3,2))/(2*mothernum+childnum*mothernum-cn(2,2));
euo=[1.7 1.4 1 1;1.7 1.5 1.3 0.3;1 1 1.4 0.3];estd=[0.1 0.1 0;0.1 0.1 0.1;0 0.1 0.1];
auo=[0 0 0 0;0 0 0 0;0 0 0 0];astd=[0.1 0.1 0;0.1 0.1 0.1;0 0.1 0.1];
for i=1:(mothernum*childnum)
if chidata(i,1)==11
        phichild(i,1:4)=[1 0 0 0];
    end
    if chidata(i,1)==12
        phichild(i,1:4)=[0 1 0 0];
    end
    if chidata(i,1)==21
        phichild(i,1:4)=[1 0 0 0];
    end
    if chidata(i,1)==22
        phichild(i,1:4)=[0 1-ep(t,1) ep(t,1) 0];
    end
    if chidata(i,1)==23
        phichild(i,1:4)=[0 0 0 1];
    end
    if chidata(i,1)==32
        phichild(i,1:4)=[0 0 1 0];
    end
    if chidata(i,1)==33
        phichild(i,1:4)=[0 0 0 1];
    end
end
box=[1 2 1 1;1 2 2 3;2 2 2 3];
while any(any(abs(auo-euo)>0.0001))|(estd(2,2)-astd)>0.0001
astd=estd;auo=euo;
for i=1:(mothernum*childnum)
    k=fix(chidata(i,1)/10);
    for j=1:4
    fpheno(i,j)=normpdf(chidata(i,2),euo(k,j),estd(fix(chidata(i,1)/10),box(fix(chidata(i,1)/10),j)));
    end
    for j=1:4
    phipheno(i,j)=phichild(i,j)*fpheno(i,j)/sum(phichild(i,:).*fpheno(i,:));
    end
end
for i=1:numel(euo(:,1))
    a=find(fix(chidata(:,1)/10)==i);
for j=1:numel(euo(1,:))
euo(i,j)=sum(chidata(a,2).*phipheno(a,j))/sum(phipheno(a,j),1);
euo(isnan(euo))=0;
minus(a,j)=phipheno(a,j).*((chidata(a,2)-euo(i,j)).^2);
end
end
for i=1:3
    a=find(fix(chidata(:,1)/10)==i);
    for j=1:3
    b=find(mod(chidata(:,1),10)==j);
    c=intersect(a,b);
    estd(i,j)=(sum(sum(minus(c,:)))/numel(c))^0.5;
    end
end
end
clear a b c astd auo ans i j k minus oricn orimn
 para_unify(t,1)=(euo(1,1)+euo(2,1)+euo(2,4)+euo(3,4))/4;
 para_unify(t,2)=(euo(1,1)+euo(2,1)-euo(2,4)-euo(3,4))/4;
 para_unify(t,3)=(sum(sum(euo))-2*(euo(1,1)+euo(2,1)+euo(2,4)+euo(3,4)))/4;
 para_unify(t,4)=(euo(1,2)+euo(2,2)-euo(2,3)-euo(3,3))/4;
 para_unify(t,5)=(euo(1,1)-euo(2,1))/2;     
 para_unify(t,6)=(euo(3,4)-euo(2,4))/2;
 para_unify(t,7)=(euo(1,2)-euo(2,2))/2;
 para_unify(t,8)=(euo(3,3)-euo(2,3))/2;
 paralogL1(t,1)=sum(log(sum(fpheno.*phichild,2)));
 clear fpheno euo phipheno
 %hypothesis test1 :all para=0%
    u0=mean(reshape(childphenotype,1,[])); std0=std(reshape(childphenotype,1,[]));
    l0=log(normpdf(reshape(childphenotype,1,[]),u0,std0));
    para1logL0(t,1)=sum(l0);
    LR(t,1)=(-2)*(para1logL0(t,1)-paralogL1(t,1));
    AIC(t,1)=(2*7-2*paralogL1(t,1))/(mothernum*childnum);
    BIC(t,1)=(7*log(mothernum*childnum)-2*paralogL1(t,1));
    clear u0 std0 l0 
  
 %hypothesis test 2: h=u1=u2=v1=v2=0 traditional model%
 tchidata(:,1)=mod(chidata(:,1),10);tchidata(:,2)=chidata(:,2);
 eu0(1)=mean(tchidata(tchidata(:,1)==1,2));
 eu0(2)=mean(tchidata(tchidata(:,1)==2,2));
 eu0(3)=mean(tchidata(tchidata(:,1)==3,2));
 estd3(1)=std(tchidata(tchidata(:,1)==1,2));
 estd3(2)=std(tchidata(tchidata(:,1)==2,2));
 estd3(3)=std(tchidata(tchidata(:,1)==3,2));
 for i=1:numel(tchidata(:,1))
     l0(i,1)=log(normpdf(tchidata(i,2),eu0(tchidata(i,1)),estd3(tchidata(i,1))));
 end
 para3logL0(t,1)=sum(l0);
 LR(t,2)=(-2)*(para3logL0(t,1)-paralogL1(t,1));
 para_tra(t,1)=0.5*(eu0(1)+eu0(3));para_tra(t,2)=0.5*(eu0(1)-eu0(3));para_tra(t,3)=0.5*(2*eu0(2)-(eu0(1)+eu0(3)));
    AIC(t,3)=(2*2-2*para3logL0(t,1))/(mothernum*childnum);
    BIC(t,3)=(2*log(mothernum*childnum)-2*para3logL0(t,1)); 
 clear eu0 estd3 l0 
 
 %hypothesis test 3: u1=u2=v1=v2=0 imprinting model%
 euo=[1 1 1 1];estd4=[0.4 0.4 0.4];
 auo=[0 0 0 0];astd=[0 0 0];
 box2=[1 2 2 3];
while any(abs(auo-euo)>0.0001)|(estd4(1,2)-astd)>0.0001
astd=estd4;auo=euo;
for i=1:(mothernum*childnum)
    for j=1:4
    fpheno(i,j)=normpdf(chidata(i,2),euo(j),estd4(box2(j)));
    end
    for j=1:4
    phipheno(i,j)=phichild(i,j)*fpheno(i,j)/sum(phichild(i,:).*fpheno(i,:));
    end
end
for j=1:4
euo(j)=sum(chidata(:,2).*phipheno(:,j))/sum(phipheno(:,j),1);
minus(:,j)=phipheno(:,j).*((chidata(:,2)-euo(1,j)).^2);
end
    for j=1:3
    b=find(mod(chidata(:,1),10)==j);
    estd4(j)=(sum(sum(minus(b,:)))/numel(b))^0.5;
    end
end
para4logL0(t,1)=sum(log(sum(fpheno.*phichild,2)));
para_imprint(t,1)=0.5*(euo(1)+euo(4));
para_imprint(t,2)=0.5*(euo(1)-euo(4));
para_imprint(t,3)=0.5*(euo(2)+euo(3)-euo(1)-euo(4));
para_imprint(t,4)=0.5*(euo(2)-euo(3));
LR(t,3)=(-2)*(para4logL0(t,1)-paralogL1(t,1));
    AIC(t,2)=(3*2-2*para4logL0(t,1))/(mothernum*childnum);
    BIC(t,2)=(3*log(mothernum*childnum)-2*para4logL0(t,1));
clear estd4 euo fpheno phipheno auo astd

%hypothesis test 4: h=0%
tchidata(:,3)=fix(chidata(:,1)/10);
for i=1:3
    for j=1:3
        euo(i,j)=mean(tchidata(tchidata(:,3)==i&tchidata(:,1)==j,2));
    end
end
 for i=1:numel(tchidata(:,1))
     l0(i,1)=log(normpdf(tchidata(i,2),euo(tchidata(i,3),tchidata(i,1)),estd(tchidata(i,3),tchidata(i,1))));
 end
 para5logL0(t,1)=sum(l0);
 LR(t,4)=(-2)*(para5logL0(t,1)-paralogL1(t,1));
 clear estd euo i j l0 minus 
 t
end

[m n]=find(isnan(paralogL1));
para_unify(m,:)=[];

Epara_unify=mean(para_unify);Epara_imprint=mean(para_imprint);Epara_tra=mean(para_tra);
Epara_unify(2,:)=std(para_unify);Epara_imprint(2,:)=std(para_imprint);Epara_tra(2,:)=std(para_tra);

for i=1:t
    AICchose(i)=find(AIC(i,:)==min(AIC(i,:)));
    BICchose(i)=find(BIC(i,:)==min(BIC(i,:)));
end
AICChose=tabulate(AICchose);
BICChose=tabulate(BICchose);

power(1)=numel(find(LR(:,1)>=14.067))/100;
power(2)=numel(find(LR(:,2)>=11.070))/100;
power(3)=numel(find(LR(:,3)>=9.488))/100;
power(4)=numel(find(LR(:,4)>=3.841))/100;
clear b box box2 chidata childgenotype childnum childphenotype cn data phichild tchidata
clear handledata i m n mn para1logL0 para3logL0 para4logL0 para5logL0 paralogL1
save('results4001unify0.2.mat');