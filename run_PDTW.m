clear all;close all;clc;

%% Data

Osignal=[]; % the oldes signal 
Odepth=[]; % Depths of oldes signal (Osignal and Odepth should have the same length) 
Ysignal=[];% the youngest signal 
Ydepth=[]; %Depths of youngest signal 

%% Settings

res=1000; % linar interpolate each signal to the length of res

w22=[Ysignal/2 Ysignal(end)]; % interval where the end of Ysignal are allowed to correspond at Osignal (m)

p=0.5; % factor P
collelengt=10; % path momentum length (relative to res)
collestrength=1; % path momentum strength
steer=0.25; % the steer factor

NN=300; % Number of simulations

constraint1=[]; % enter constraints. E.g if a depth within 50-60 meter at site1 correspond to site2 60-70 m then enter constraint1=[50 60 60 70], add another constraint by e.g. constraint1=[50 60 60 70; 1 2 3 4], etc.

nd=0; % if nd=1 then plot the dtw alignment


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% code start from here

[~,w3]=min(abs(Odepth'-w22));w2=length(Odepth)-w3;


if length(constraint1)>0
constraint=zeros(size(constraint1));
for i=1:size(constraint,1);
[~,constraint(i,1:2)]=min(abs(Odepth(:)-constraint1(i,1:2)));[~,constraint(i,3:4)]=min(abs(Ysingal(:)-constraint1(i,3:4)));end
else
constraint=0;
end

endB=1; % if 0 then end at ns,nt
w=1; % define cutoff on edges

L1=Odepth(end);L2=Ydepth(end); % interpolation
a2=interp1(Odepth(:).*linspace(1,1.00000000001,length(Odepth))',Osignal,linspace(0,L1,res))';a=a2(a2>0);Odepth=linspace(0,L1,res);Odepth=Odepth(a2>0); 
b2=interp1(Ydepth(:).*linspace(1,1.00000000001,length(Ydepth))',Ysignal,linspace(0,L2,res))';b=b2(b2>0);Ydepth=linspace(0,L2,res);Ydepth=Ydepth(b2>0);

[d1 DD paths k]=PDTW(a,b,w,w2,p,NN,collelengt,collestrength,steer,constraint);


%% plots
% dencity plot

h(1)=subplot(221);
pa2=[];pb2=[];
for i=1:NN
    pa1=paths(paths(:,1,i)>1 ,1,i); %& paths(:,1,i)<length(a) & paths(:,2,i)<length(b)
    pb1=paths(paths(:,2,i)>1 ,2,i); %& paths(:,1,i)<length(a)
    alignments{i}=[Ydepth(pb2-1);Odepth(pa2-1)];
    pa2=[pa2;pa1];pb2=[pb2;pb1];%
end
    NX=201;NY=201;
    [C,x_arr1,y_arr1]=histcounts2(Ydepth(pb2-1),Odepth(pa2-1),[NX,NY]);
    x_arr=(x_arr1(1:end-1)+x_arr1(2:end))/2;y_arr=(y_arr1(1:end-1)+y_arr1(2:end))/2;
    P=C./sum(C(:)); % THIS IS THE 2D PRIOR PDF!!
    imagesc(x_arr,y_arr,P');caxis([0 0.001]);colormap(1-gray); colorbar; box on;grid on;

    if nd==1
[d12 DD2 thepath k2]=PDTW(a,b,w,w2,1,1,0,0,0,0);%(a,b,w,w2,p,NN,collelengt,collestrength,styr,likelihood,constraint);
pa=thepath(thepath(:,1,1)>0 & thepath(:,1,1)<length(a) & thepath(:,2,1)<length(b),1,1); pb=thepath(thepath(:,2,1)>0 & thepath(:,2,1)<length(b) & thepath(:,1,1)<length(a),2,1);
hold on;plot(Ydepth(pb),Odepth(pa),'b',LineWidth=1); 
    end

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

PDTW_pdf=P;
save('P-DTW_results.mat','alignments','PDTW_pdf');

