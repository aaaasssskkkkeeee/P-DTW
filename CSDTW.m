
function [d D paths k]=dtw(a,b,w,w2,p,NN,collelengt,collestrength,styr,likelihood,constraint)


if nargin<3
    w=Inf;
end
na=size(a,1); % nr of data points in a
nb=size(b,1);

w=max(w,abs(na-nb)); % window size


%% initialization

numbers(1)=3;
paths=zeros((na+1)+(nb+1),2,NN);

%% Calculate the cost matrix

if likelihood==1
D=zeros(na+1,nb+1); 
D(1,1)=1;
errora=a*p;errorb=a*p;
    for i=1:na
        for j=max(i-w,1):min(i+w,nb)
            logL=exp(-0.5*(a(i,:)-b(j,:))^2/(errora(i)^2+errorb(j)^2));
            D(i+1,j+1)=logL*max( [D(i,j+1), D(i+1,j), D(i,j)] ); 
        end
    end
else

D=zeros(na+1,nb+1); 
D(:,1)=Inf;D(1,:)=Inf;
D(1,1)=0;

if length(constraint)>3 % add constraints
    for i=1:size(constraint,1)
        D(1:constraint(i,1)+1,constraint(i,4)+1:end)=inf;
        D(constraint(i,2)+1:end,1:constraint(i,3)+1)=inf;
    end
end
    for i=1:na % calculate cost matrix
        for j=1:nb
            if D(i+1,j+1)<inf
            resi=abs(a(i)-b(j)); %* (errora(i)^2+errorb(j)^2)^0.5;
            D(i+1,j+1)=resi + min( [D(i,j+1), D(i+1,j), D(i,j)] );
            end
        end
    end
end



for iNN=1:NN
%% Randomly choose initial starting point

if likelihood==1
    %ag=[D((ns-w2+2):ns+1,nt+1)./((ns-w2+2):ns+1)'; D(ns+1,(nt-w2+2):nt+1)'./((nt-w2+2):nt+1)'];
    ag=[D((na-w2(1)+2):na+1-w2(2),nb+1)./((na-w2(1)+2):(na+1-w2(2)))'];
    M=ag.^0.5; %-ag(nr(1));
else

    agI=[D((na-w2(1)+2):na+1-w2(2),nb+1)];%./((na-w2(1)+2):na+1-w2(2))';
    [~, nrI]=sort(agI);
    M=(1./(agI - agI(nrI(1)) + (median(agI)-agI(nrI(1)))*(1-p)));

end


index=1+sum([cumsum(M)]/sum(M)<rand);  
d=M(index); 

if index<=w2(1)-w2(2)
    n=(na-w2(1))+index;
    m=nb+1;
else
    n=na+1;
    m=(nb-w2+1)+(index-w2);
end


%% find path

path=[];
k=1;
path(1,:)=[n,m];
numbers(1)=1;

while ((n+m)~=2) % stops when n=1 & m=1
    
    if (n-1)==0
        m=m-1;
    elseif (m-1)==0
        n=n-1;
    else 

            ag=[D(n-1,m),D(n,m-1),D(n-1,m-1)];
   
   
   [~, nr]=sort(ag);
   if ag(nr(2)) == inf
       number = nr(1);
   else

           styring=[n/m m/n 1].^styr; % go towards end if all 3 steps are eqally probable
   
        if likelihood==1
            M=ag.*styring;
        else
            M=(1./(ag - ag(nr(1))+(ag(nr(2))-ag(nr(1)))*(1-p))).*styring; % relative to secund best
        end
       
       
        if collestrength>0
            colst=[(path(max(k-collelengt,1),1)-n+1)/(path(max(k-collelengt,1),2)-m+1) (path(max(k-collelengt,1),2)-m+1)/(path(max(k-collelengt,1),1)-n+1) 1];
            M=M.*colst.^collestrength;
        end

       number=1+sum([M(1), M(1)+M(2)]/sum(M)<rand); % choose path based on probability
  end

      switch number
      case 1
        n=n-1;
      case 2
        m=m-1;
      case 3
        n=n-1;
        m=m-1;
      end
  end
    k=k+1;
    path=cat(1,path,[n,m]);
    numbers=cat(1,numbers,number);
end
clear numbers

paths(1:size(path,1),:,iNN)=path;

end

