clear all;close all;clc;

%% Data

Osignal=[]; % the oldes signal 
Odepth=[];% Depth of oldes signal 
Ysingal=[];% the youngest signal 
Ydepth=[]; %Depth of youngest signal 

%% Settings

res=4000; % linar interpolate each signal to the length of res

    w22=[50 100];
    [~,w3]=min(abs(Odepth'-w22));w2=length(Odepth)-w3;


p=0.5; % factor P
collelengt=10; % path momentum length
collestrength=1; % path momentum strength
steer=0.25; % the steer factor

NN=100; % Number of alignments

constraint1=[]; % enter constraints. E.g if a depth within 50-60 meter at site1 correspond to site2 60-70 m then enter constraint1=[50 60 60 70], add another constraint by e.g. constraint1=[50 60 60 70; 1 2 3 4], etc.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% code start from here
if length(constraint1)>0
constraint=zeros(size(constraint1));
for i=1:size(constraint,1);
[~,constraint(i,1:2)]=min(abs(depthd181'-constraint1(i,1:2)));[~,constraint(i,3:4)]=min(abs(depthd182'-constraint1(i,3:4)));end
else
constraint=0;
end

endB=1; % if 0 then end at ns,nt
w=1; % define cutoff on edges
likelihood=0;

L1=Odepth(end);L2=Ydepth(end); % interpolation
a2=interp1(Odepth.*linspace(1,1.00000000001,length(Odepth))',a,linspace(0,L1,res))';a=a2(a2>0);Odepth=linspace(0,L1,res);Odepth=Odepth(a2>0); 
b2=interp1(Ydepth.*linspace(1,1.00000000001,length(Ydepth))',b,linspace(0,L2,res))';b=b2(b2>0);Ydepth=linspace(0,L2,res);Ydepth=Ydepth(b2>0);

[d1 DD paths k]=CSDTW(a,b,w,w2,p,NN,collelengt,collestrength,steer,likelihood,constraint);


%% plots
% dencity plot

h(1)=subplot(221);
pa2=[];pb2=[];
for i=1:NN
    pa1=paths(paths(:,1,i)>1 ,1,i); %& paths(:,1,i)<length(a) & paths(:,2,i)<length(b)
    pb1=paths(paths(:,2,i)>1 ,2,i); %& paths(:,1,i)<length(a)
    pa2=[pa2;pa1];pb2=[pb2;pb1];%
end
    NX=201;NY=201;
    [C,x_arr1,y_arr1]=histcounts2(depthd182(pb2-1),depthd181(pa2-1),[NX,NY]);
    x_arr=(x_arr1(1:end-1)+x_arr1(2:end))/2;y_arr=(y_arr1(1:end-1)+y_arr1(2:end))/2;
    P=C./sum(C(:)); % THIS IS THE 2D PRIOR PDF!!
    imagesc(x_arr,y_arr,P');caxis([0 0.001]);colormap(1-gray); colorbar; box on;grid on;

% det=interp1(Aage2(depthd182),depthd182,Aage1(depthd181)); % real depth2 as a f og depth1
% hold on;plot(det(det>0),depthd181(det>0),'g','LineWidth',1)

[d12 DD2 thepath k2]=sdtw(a,b,w,w2,1,1,0,0,0,0,0);%(a,b,w,w2,p,NN,collelengt,collestrength,styr,likelihood,constraint);
pa=thepath(thepath(:,1,1)>0 & thepath(:,1,1)<length(a) & thepath(:,2,1)<length(b),1,1); pb=thepath(thepath(:,2,1)>0 & thepath(:,2,1)<length(b) & thepath(:,1,1)<length(a),2,1);
hold on;plot(Ydepth(pb),Odepth(pa),'b',LineWidth=1); 
if constraint>0;hold on;plot(constraint1(:,4:-1:3),constraint1(:,1:2),'*r');end

set(gca,'Ydir','reverse','Xdir','reverse');
ylim([0 Odepth(end)]);xlim([0 Ydepth(end)]);
grid on; box on;

h(2)=subplot(222);plot(a, Odepth,'k');ylabel('Depth (m)');box on; grid on;

title('site 1');
set(gca,'Ydir','reverse');ylim([0 Odepth(end)]);

h(3)=subplot(223);plot(Ydepth,b,'k');set(gca,'Xdir','reverse');xlabel('Depth (m)');box on; grid on;
title('site 2');
xlim([0 Ydepth(end)]);
set(h(1),'Position',[0.3-0.1  0.3-0.1 0.6+0.1  0.6+0.1 ],'fontsize',10,'xticklabel',[],'yticklabel',[]);set(h(2),'Position',[0.1    0.3-0.1    0.2-0.1    0.6+0.1 ],'fontsize',10);set(h(3),'Position',[0.3-0.1     0.1    0.6+0.1    0.2-0.1 ],'fontsize',10);


%% line correlation plot
places=[10:20:90]; % depths at the youngest site to correlate

ost=zeros(length(places),NN);
y=[];x=[];
for i=1:NN
    pa1=paths(paths(:,1,i)>1 ,1,i); %& paths(:,1,i)<length(a) & paths(:,2,i)<length(b)
    pb1=paths(paths(:,2,i)>1 ,2,i); %& paths(:,1,i)<length(a)
     der=interp1(depthd181(pa1-1)+linspace(0,1e-10,length(pb1)),depthd182(pb1-1)+linspace(0,1e-10,length(pb1)),places);
    ost(:,i)=der;
end
    
    y=sort(rand(1000,NN));
    X=[];Y=[];
for i=1:NN
    x= (places'-ost(:,i)).*y(:,i)'+places';
    X=[X x(:)];
    yy=repmat(y(:,i)',length(places),1);
    Y=[Y yy(:)];
end
NX=201;NY=201;
    [C,x_arr1,y_arr1]=histcounts2(X(:),Y(:),[NX,NY]);
    x_arr=(x_arr1(1:end-1)+x_arr1(2:end))/2;y_arr=(y_arr1(1:end-1)+y_arr1(2:end))/2;
    P=C./sum(C(:)); 

h(1)=subplot(311);plot(depthd182,b,'k');box on; grid on;xlim([0 depthd181(end)])
%hold on;xline(places)
h(2)=subplot(312);
    imagesc(x_arr,y_arr,P');caxis([0 0.001]);colormap(1-gray); 
    xlim([0 depthd181(end)]);box off
h(3)=subplot(313);
plot(depthd181,a,'k');xlabel('Depth (m)');box on; grid on;xlim([0 depthd181(end)])
set(h(1),'Position',[0.1  0.7 0.8  0.25],'fontsize',10,'xticklabel',[]);
set(h(2),'Position',[0.1   0.7-0.25    0.8   0.25 ],'fontsize',10,'xticklabel',[],'ytick',[],'YColor','none');
set(h(3),'Position',[0.1     0.2    0.8    0.25 ],'fontsize',10);
