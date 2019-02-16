clc
close all
clear all
ds = datastore('house_data_complete.csv','TreatAsMissing','NA',.....
    'MissingValue',0,'ReadSize',25000);
T = read(ds);
size(T);
Alpha=.005;
m=length(T{:,1});
p=ceil(0.6*m);
p1=ceil(0.2*m);
U0=T{1:p ,2}
U=T{1:p,4:19};
U1=T{1:p,20:21};
U2=U.^2;
U3=U.^4;
X=[ones(p,1) U U1 U3]

n=length(X(1,:));
for w=2:n
    if max(abs(X(:,w)))~=0
    X(:,w)=(X(:,w)-mean((X(:,w))))./std(X(:,w));
    end
end

Y=T{1:p,3}/mean(T{1:p,3});
Theta=zeros(n,1);

k=1;

E(k)=(1/(2*(p)))*sum((X*Theta-Y).^2);

R=1;
while R==1
Alpha=Alpha*1;
Theta=Theta-(Alpha/p)*X'*(X*Theta-Y);
k=k+1
E(k)=(1/(2*(p))*sum((X*Theta-Y).^2));
if E(k-1)-E(k)<0
    break
end 
q=(E(k-1)-E(k))./E(k-1);
if q <.0001;
    R=0;
end
end
plot(E);

%%%%%%%


U00=T{p+1:p1+p ,2}
UU=T{p+1:p1+p,4:19};
U11=T{p+1:p1+p,20:21};
U22=UU.^2;
U33=UU.^4;
X1=[ones(p1,1) UU U11 U33]
Theta1=Theta;
YY=T{p+1:p1+p,3}/mean(T{p+1:p1+p,3});

k=1;
n1=length(X1(1,:));
for w=2:n1
    if max(abs(X1(:,w)))~=0
    X1(:,w)=(X1(:,w)-mean((X1(:,w))))./std(X1(:,w));
    end
end
E1(k)=(1/(2*(p1)))*sum((X1*Theta1-YY).^2);


%%%%%%%

U000=T{p+p1+1:end ,2}
UUU=T{p+p1+1:end,4:19};
U111=T{p+p1+1:end,20:21};
U222=UUU.^2;
U333=UUU.^4;
X2=[ones(length(UUU),1) UUU U111 U333]
Theta2=Theta1;
YYY=T{p+p1+1:end,3}/mean(T{p+p1+1:end,3});

k=1;
n2=length(X2(1,:));
for w=2:n2
    if max(abs(X2(:,w)))~=0
    X2(:,w)=(X2(:,w)-mean((X2(:,w))))./std(X2(:,w));
    end
end
E2(k)=(1/(2*(p1)))*sum((X2*Theta2-YYY).^2);




