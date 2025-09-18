IGD_m=[
    0.1204    0.4058    0.4253    0.5197    0.4816
    0.2600    0.4054    0.5876    0.5293    0.6411
    0.1250    0.2659    0.5765    0.5197    0.7286
    0.1533    0.3296    0.1024    0.5302    0.4968
    0.2571    0.4351    0.4256    0.4982    0.4982
    0.1381    0.4163    0.4791    0.3291    0.6300
    0.1489    0.2359    0.6237    0.5509    0.6652
    0.1486    0.2163    0.5305    0.5497    0.5787
    0.1512    0.2177    0.4518    0.5205    0.8495
    0.1372    0.4276    0.2113    0.4982    0.7284
    0.1205    0.2781    0.6876    0.5565    0.6402
    0.2694    0.4403    0.4272    0.5151    0.7298
    0.2732    0.3221    0.5220    0.5113    0.7295
    0.2962    0.4598    0.6151    0.5197    0.4194
    0.1720    0.4079    0.7824    0.5797    0.6654
    0.0937    0.6776    0.6073    0.5293    0.5215
    0.2112    0.3281    0.5016    0.5197    0.6300
    0.2739    0.3266    0.5207    0.5197    0.6654
    0.1856    0.2885    0.6847    0.5409    0.7286
    0.1361    0.3087    0.8655    0.5675    0.9559
];

HV_m=[
    0.4337    0.2619    0.3538    0.3387    0.0158
    0.6463    0.2998    0.1815    0.4079    0.0517
    0.5594    0.3388    0.2982    0.4281    0.0458
    0.4494    0.2331    0.1352    0.4063    0.0271
    0.4919    0.2838    0.0647    0.4238    0.0795
    0.4721    0.2482    0.2485    0.2617    0.0538
    0.5316    0.2329    0.2442    0.2515    0.0475
    0.4924    0.0915    0.2867    0.0600    0.0221
    0.5957    0.1454    0.2306    0.4536    0.0243
    0.5982    0.3100    0.2431    0.1712    0.0459
    0.5929    0.2307    0.0901    0.4162    0.0518
    0.5729    0.2467    0.4412    0.3065    0.0208
    0.5312    0.2848    0.0128    0.4571    0.0457
    0.5959    0.3986    0.4215    0.4452    0.2594
    0.4481    0.0481    0.2426    0.4390    0.0475
    0.5037    0.0841    0.3730    0.4014    0.0910
    0.5338    0.3217    0.1554    0.2420    0.0538
    0.6115    0.2167    0.3750    0.4547    0.0474
    0.6106    0.1091    0.2658    0.3397    0.0458
    0.6025    0.3532    0.3151    0.3564    0.0226];
Fddd=min(IGD_m);     %IGD最小值
Fddd1=mean(IGD_m);   %IGD平均值
Feee=max(HV_m);      %HV最大值
Feee1=mean(HV_m);    %HV平均值

zhibiao111=[Fddd;Fddd1;Feee;Feee1];

HV_3=HV_m;

[~,a]=min(min(IGD_m));
[~,a1]=min(mean(IGD_m));

[~,c]=max(max(HV_m));
[~,c1]=max(mean(HV_m));

IGD_m   
HV_m
%% 处理IGD
means_igd = mean(IGD_m )
mins_igd = min(IGD_m )


%% 处理HV
means_hv = mean(HV_m)
max_hv = max(HV_m)


[a,a1,c,c1]    %排名

[~,ww1]=sort(min(IGD_m));
[~,ww2]=sort(mean(IGD_m));
[~,ww3]=sort(max(HV_m));
[~,ww4]=sort(mean(HV_m));

format long   %kstest2
[h1,p1,k]=ranksum(HV_m(:,1),HV_m(:,2),0.05);
[h2,p2,k]=ranksum(HV_m(:,1),HV_m(:,3),0.05);
[h3,p3,k]=ranksum(HV_m(:,1),HV_m(:,4),0.05);
[h4,p4,k]=ranksum(HV_m(:,1),HV_m(:,5),0.05);
w9=[p1,p2,p3,p4];
m9=[h1,h2,h3,h4];
kHV=double([m9])

[h1,p1,k]=ranksum(IGD_m(:,1),IGD_m(:,2),0.05);
[h2,p2,k]=ranksum(IGD_m(:,1),IGD_m(:,3),0.05);
[h3,p3,k]=ranksum(IGD_m(:,1),IGD_m(:,4),0.05);
[h4,p4,k]=ranksum(IGD_m(:,1),IGD_m(:,5),0.05);
w9=[p1,p2,p3,p4];
z9=[h1,h2,h3,h4];
kIGD=double([z9])

zhibiao_last1(:,1)=IGD_m(:,1);
zhibiao_last1(:,3)=IGD_m(:,2);
zhibiao_last1(:,5)=IGD_m(:,3);
zhibiao_last1(:,7)=IGD_m(:,4);
zhibiao_last1(:,9)=IGD_m(:,5);

zhibiao_last1(:,2)=HV_m(:,1);
zhibiao_last1(:,4)=HV_m(:,2);
zhibiao_last1(:,6)=HV_m(:,3);
zhibiao_last1(:,8)=HV_m(:,4);
zhibiao_last1(:,10)=HV_m(:,5);

figure 
zhibiao=zhibiao_last1;
IGD=zeros(20,5);
IGD(:,1)=zhibiao(:,1);
IGD(:,2)=zhibiao(:,3);
IGD(:,3)=zhibiao(:,5);
IGD(:,4)=zhibiao(:,7);
IGD(:,5)=zhibiao(:,9);

HV=zeros(20,5);
HV(:,1)=zhibiao(:,2);
HV(:,2)=zhibiao(:,4);
HV(:,3)=zhibiao(:,6);
HV(:,4)=zhibiao(:,8);
HV(:,5)=zhibiao(:,10);

subplot(1,2,1)
h=boxplot(IGD);
g = gca;
line=g.Children.Children;
line(1).Color ='b';
line(2).Color ='b';
line(3).Color ='b';
line(4).Color ='b'; %最下面的线
line(5).Color ='b';
line(6).Color ='r';
line(7).Color ='r';
line(8).Color ='r';
line(9).Color ='r';
line(10).Color ='r';
line(11).Color ='b'; %最下面的线
line(12).Color ='b';
line(13).Color ='b';
line(14).Color ='b';
line(15).Color ='b';
line(16).Color ='k';
line(17).Color ='k';
line(18).Color ='k'; %最下面的线
line(19).Color ='k';
line(20).Color ='k';
line(21).Color ='k';
set(h,'Linewidth',1.5');
set(gca,'XTickLabel',{'QDS-KOA','KOA','NSGAII','VNS-NSGAII','GSA'})
ylabel('IGD')

subplot(1,2,2)
n=boxplot(HV);
g2 = gca;
line1=g2.Children.Children;
line(1).Color ='b';
line(2).Color ='b';
line(3).Color ='b';
line(4).Color ='b'; %最下面的线
line(5).Color ='b';
line(6).Color ='r';
line(7).Color ='r';
line(8).Color ='r';
line(9).Color ='r';
line(10).Color ='r';
line(11).Color ='b'; %最下面的线
line(12).Color ='b';
line(13).Color ='b';
line(14).Color ='b';
line(15).Color ='b';
line(16).Color ='k';
line(17).Color ='k';
line(18).Color ='k'; %最下面的线
line(19).Color ='k';
line(20).Color ='k';
line(21).Color ='k';
set(n,'Linewidth',1.5');
set(gca,'XTickLabel',{'QDS-KOA','KOA','NSGAII','VNS-NSGAII','GSA'})
% set(gca,'YTickLabel','S');
ylabel('HV')


