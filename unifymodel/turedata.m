    clc
    clear
    cd ('F:\MATLAB\work\');
    motherdata=xlsread('torreya_srap_data.xls','parent_data');
    childdata=xlsread('torreya_srap_data.xls','child_data');
    phenotypedata=xlsread('torreya_srap_data.xls','phenotype');
    a=childdata(:,3:235);a(isnan(a))=1;a(a~=0)=1;
    childdata(:,3:235)=a;
    clear a
    childnum=20;jj=3;
    for ii=3:4
         %handle parent_data%
   phenodata=phenotypedata(:,jj);
   momdata=motherdata(:,ii);
   chidata=childdata(:,ii);
   for i=1:numel(momdata)
       chidata(childnum*i-19:childnum*i,2)=momdata(i,1);
   end
   mn(1,1)=numel(find(momdata==1));  mn(1,2)=numel(find(momdata==0));
   cn(1,1)=numel(find(chidata(:,1)==1&chidata(:,2)==1));
   cn(1,2)=numel(find(chidata(:,1)==0&chidata(:,2)==1));
   cn(2,1)=numel(find(chidata(:,1)==1&chidata(:,2)==0));
   cn(2,2)=numel(find(chidata(:,1)==0&chidata(:,2)==0));
   
   epm=0.8;epo=0.9;apm=0;apo=0;phi=[0.7 0.6 0.5];
   while abs(epm-apm)>0.0001|abs(epo-apo)>0.0001
   apm=epm;apo=epo;
   epm=(phi(1)*mn(1,1)+phi(2)*cn(1,1))/sum(2*sum(mn)+sum(sum(cn)));
   epo=(cn(2,1)+phi(3)*cn(1,1))/sum(sum(cn));
   phi(1)=(2*epm^2+2*epm*(1-epm))/(epm^2+2*epm*(1-epm));
   phi(2)=(epm^2+epm*(1-epm))/(epm^2+epm*(1-epm)*(1+epo));
   phi(3)=(epm^2*epo+2*epm*(1-epm)*epo)/(epm^2+epm*(1-epm)*(1+epo));
   end
   
   Ep(ii,1)=epm;Ep(ii,2)=epo;
    for i=1:numel(momdata)
       if any(chidata(20*i-19:20*i)==0)&momdata(i,1)==1
           momdata(i,1)=2;
       end
       chidata(childnum*i-19:childnum*i,2)=momdata(i,1);
   end
   orimn=tabulate(momdata(:));
   mn=zeros(3,1);
   for i=1:3
       if isempty(orimn(orimn(:,1)==(i-1),2))
           continue;
       end
       mn(i,1)=orimn(orimn(:,1)==(i-1),2);
   end
   mothernum=sum(mn);
    %estimate L0%
    chidata(isnan(phenodata),:)=[];phenodata(isnan(phenodata))=[];
    u0=mean(phenodata); std0=std(phenodata);
    l0=log(normpdf(phenodata,u0,std0));
    L0(jj,ii)=sum(l0);
    
     %L1 for traditional model%
       euo=[u0 u0+1 u0-1];estd=std0;
       auo=[0 0 0];astd=0;
       phichild=zeros(numel(chidata(:,1)),3);      
       phichild(chidata(:,1)==0,1)=1;
phichild(chidata(:,1)==1,2)=(epo*epm^2+epo*epm*(1-epm))/(epm^2+epm*(1-epm)*(1+epo));
phichild(chidata(:,1)==1,3)=1-(epo*epm^2+epo*epm*(1-epm))/(epm^2+epm*(1-epm)*(1+epo));
while any(abs(auo-euo)>0.001)|abs(estd-astd)>0.001
astd=estd;auo=euo;
    for j=1:3
    fpheno(:,j)=normpdf(phenodata(:,1),euo(j),estd);
    end
    for j=1:3
    phipheno(:,j)=phichild(:,j).*fpheno(:,j)./sum(phichild.*fpheno,2);
    euo(j)=sum(phenodata(:,1).*phipheno(:,j))/sum(phipheno(:,j),1);
    minus(:,j)=phenodata(:,1)-euo(j);
    end
estd=(sum(sum((minus.^2).*phipheno))/numel(chidata(:,1)))^0.5;
end
tL1(jj,ii)=sum(log(sum(fpheno.*phichild,2)));
 para_tra(ii,1)=0.5*(euo(1)+euo(2));para_tra(ii,2)=0.5*(euo(2)-euo(1));
 para_tra(ii,3)=0.5*(2*euo(3)-(euo(1)+euo(2)));
 AIC(ii,3)=(2*2-2*L1(jj,ii))/(numel(phenodata));
clear auo astd phichild phipheno

 %L1 for impringting model%
    iuo=[u0 u0+1+0.5 u0+1-0.5 u0-1];istd=std0;
    auo=[0 0 0 0];astd=0;
    phichild=zeros(numel(chidata(:,1)),4);   
    phichild(chidata(:,1)==0&chidata(:,2)==0,4)=1;
phichild(chidata(:,1)==1&chidata(:,2)==0,3)=1;
phichild(chidata(:,1)==0&chidata(:,2)==2,4)=1;
phichild(chidata(:,1)==1&chidata(:,2)==2,1)=epo/(2*epo+(1-epo));
phichild(chidata(:,1)==1&chidata(:,2)==2,2)=(1-epo)/(2*epo+(1-epo));
phichild(chidata(:,1)==1&chidata(:,2)==2,3)=epo/(2*epo+(1-epo));
phichild(chidata(:,1)==1&chidata(:,2)==1,1)=epo;
phichild(chidata(:,1)==1&chidata(:,2)==1,2)=1-epo;
while any(abs(auo-iuo)>0.001)|abs((istd-astd))>0.001
astd=istd;auo=iuo;
    for j=1:4
    fpheno(:,j)=normpdf(phenodata(:,1),iuo(j),istd);
    end
    for j=1:4
    phipheno(:,j)=phichild(:,j).*fpheno(:,j)./sum(phichild.*fpheno,2);
    iuo(j)=sum(phenodata(:,1).*phipheno(:,j))/sum(phipheno(:,j),1);
    minus(:,j)=phenodata(:,1)-iuo(j);
    end
istd=(sum(sum((minus.^2).*phipheno))/numel(chidata(:,1)))^0.5;
end
para_imprint(ii,1)=0.5*(iuo(1)+iuo(4));
para_imprint(ii,2)=0.5*(iuo(1)-iuo(4));
para_imprint(ii,3)=0.5*(iuo(2)+iuo(3)-iuo(1)-iuo(4));
para_imprint(ii,4)=0.5*(iuo(2)-iuo(3));
iL1(jj,ii)=sum(log(sum(fpheno.*phichild,2)));
AIC(ii,2)=(3*2-2*iL1(jj,ii))/(numel(phenodata));
  clear auo astd phichild phipheno
  
   %L1 for unify model%
    ieuo=[u0 u0+1+0.5 u0+1-0.5 u0-1;u0 u0+1+0.5 u0+1-0.5 u0-1;u0 u0+1+0.5 u0+1-0.5 u0-1];iestd=[std0 std0;std0 std0;std0 std0];
    iauo=[0 0 0 0;0 0 0 0;0 0 0 0];iastd=[0.1 0;0.1 0.1;0.1 0.1];
    box=[1 1 1 1;1 1 1 2;1 1 1 2];
    chidata(chidata(:,2)==0,2)=3; 
    chidata(chidata(:,1)==0,1)=10*chidata(chidata(:,1)==0,2)+2;
    chidata(chidata(:,1)==1,1)=10*chidata(chidata(:,1)==1,2)+1;
    chidata(:,2)=phenodata;
      phichild=zeros(numel(chidata(:,1)),4);    
for i=1:numel(chidata(:,1))
    if chidata(i,1)==11
        phichild(i,1:4)=[epo 1-epo 0 0];
    end
    if chidata(i,1)==21
        phichild(i,1:4)=[epo/(1+epo) (1-epo)/(1+epo) epo/(1+epo) 0];
    end
    if chidata(i,1)==22
        phichild(i,1:4)=[0 0 0 1];
    end
    if chidata(i,1)==31
        phichild(i,1:4)=[0 0 1 0];
    end
    if chidata(i,1)==32
        phichild(i,1:4)=[0 0 0 1];
    end
end
    
while any(any(abs(iauo-ieuo)>0.0001))|(iestd(2,2)-iastd)>0.0001
iastd=iestd;iauo=ieuo;
for i=1:numel(chidata(:,1))
    k=fix(chidata(i,1)/10);
    for j=1:4
      fpheno(i,j)=normpdf(chidata(i,2),ieuo(k,j),iestd(fix(chidata(i,1)/10),box(fix(chidata(i,1)/10),j)));
    end
    for j=1:4
    phipheno(i,j)=phichild(i,j)*fpheno(i,j)/sum(phichild(i,:).*fpheno(i,:));
    end
end
for i=1:numel(ieuo(:,1))
    a=find(fix(chidata(:,1)/10)==i);
for j=1:numel(ieuo(1,:))
ieuo(i,j)=sum(chidata(a,2).*phipheno(a,j))/sum(phipheno(a,j),1);
ieuo(isnan(ieuo))=0;
minus(a,j)=phipheno(a,j).*((chidata(a,2)-ieuo(i,j)).^2);
end
end
for i=1:3
    a=find(fix(chidata(:,1)/10)==i);
    for j=1:3
    b=find(mod(chidata(:,1),10)==j);
    c=intersect(a,b);
    iestd(i,j)=(sum(sum(minus(c,:)))/numel(c))^0.5;
    end
end
end
clear a b c astd auo 
 para_unify(ii,1)=(ieuo(1,1)+ieuo(2,1)+ieuo(2,4)+ieuo(3,4))/4;
 para_unify(ii,2)=(ieuo(1,1)+ieuo(2,1)-ieuo(2,4)-ieuo(3,4))/4;
 para_unify(ii,3)=(sum(sum(ieuo))-2*(ieuo(1,1)+ieuo(2,1)+ieuo(2,4)+ieuo(3,4)))/4;
 para_unify(ii,4)=(ieuo(1,2)+ieuo(2,2)-ieuo(2,3)-ieuo(3,3))/4;
 para_unify(ii,5)=(ieuo(1,1)-ieuo(2,1))/2;     
 para_unify(ii,6)=(ieuo(3,4)-ieuo(2,4))/2;
 para_unify(ii,7)=(ieuo(1,2)-ieuo(2,2))/2;
 para_unify(ii,8)=(ieuo(3,3)-ieuo(2,3))/2;
 uL1(jj,ii)=sum(log(sum(fpheno.*phichild,2)));
 AIC(ii,1)=(2*7-2*uL1(jj,ii))/(numel(phenodata));
 LR(ii,1)=(-2)*(L0(jj,ii)-tL1(jj,ii));
 LR(ii,2)=(-2)*(L0(jj,ii)-iL1(jj,ii));
 LR(ii,3)=(-2)*(L0(jj,ii)-uL1(jj,ii));
 ii
  clear astd auo ep eq i j mn minus momdata orimn l0 estd euo fpheno ans  u0 std0 phichild phipheno phenodata
  clear iastd iauo iestd ieuo chidata apm apo epm epo istd iuo phi 
    end
    
 clear ii jj k cn box
    
save('tdata_analysis_d2.mat')