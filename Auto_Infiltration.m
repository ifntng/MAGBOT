%% In this code, the following script parameters are required:
%WindowAPI This script is used to display images at a specified position on the screen.
%https://ww2.mathworks.cn/matlabcentral/fileexchange/31437-windowapi
%cameraParams_MAGBOT.mat This parameter is used for pre-calibrating the camera globally.
if exist('vid','var')
    delete(vid);
end

close all;
clearvars -except BG;
warning('off');

global cameraParams MAINFIG
load cameraParams_MAGBOT.mat; 

FigH   = figure(1);
axes('Visible', 'off', 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
WindowAPI(FigH, 'Position', [2405, 750, 360,360], 2); 
WindowAPI(FigH, 'Clip');

STphoto=imread('Outer frame base map.png');

P=ones(360,360,3);
gl2=0.1;
P=STphoto*gl2;

MAINFIG=image(P);
drawnow;

pphoto=cell(1,1);

MAINFIG=image(P);
drawnow

vid =videoinput('winvideo', 1, 'MJPG_1920x1080');
set(vid,'TriggerRepeat',Inf);
set(vid,'FramesPerTrigger',1);
start(vid);

for step=1:10000
    tic
    
    step
    gl1=1;
    LSD=0.5;
    P=STphoto*gl2;    
    
    P_ini=getsnapshot(vid);
    P_cut=P_ini(15:1065,460:1510,:);
    P_rgb=rgb2gray(P_cut);
    P_TW = im2bw(P_rgb,0.7);
    se = strel('disk',3);
    P_close = imclose(P_TW,se);
    
    flushdata(vid);
    
    stats=regionprops(P_close,'Area','BoundingBox','Centroid', 'MajorAxisLength'); 
    
    centers=[];
    for b=1:size(stats,1)
            centers(b,1)=round(getfield(stats,{b,1},'Centroid',{1}));
            centers(b,2)=round(getfield(stats,{b,1},'Centroid',{2})); 
    end
    
    pphoto{step,1}=centers;
    
    centers(:,1)=round((centers(:,1))*0.36-7);
    centers(:,2)=round((centers(:,2))*0.36+10);
    centersxy=[];

    centersX=round(500*0.36+2);
    centersY=round(500*0.36+2);
    
    centersxy(:,1)=(centers(:,1)-centersX);
    centersxy(:,2)=(centers(:,2)-centersY);
    
    for m=1:size(centersxy,1)
        LX(m,1)=centersxy(m,1)/centersX;
        LY(m,1)=centersxy(m,2)/centersY; 
        if centers(m,1)<centersX && centers(m,2)<centersY 
            centers(m,1)=centers(m,1)-round(LX(m,1)*10);
            centers(m,2)=centers(m,2)-round(LY(m,1)*10)-3;
            continue;
        end
        if centers(m,1)>=centersX && centers(m,2)<centersY
            centers(m,1)=centers(m,1)-round(LX(m,1)*10)-1;
            centers(m,2)=centers(m,2)-round(LY(m,1)*10-3)-8;
            continue;
        end
        if centers(m,1)<centersX && centers(m,2)>=centersY
            centers(m,1)=centers(m,1)-round(LX(m,1)*10)+2;
            centers(m,2)=centers(m,2)-round(LY(m,1)*10)-5;
            continue;
        end
        if centers(m,1)>=centersX && centers(m,2)>=centersY
            centers(m,1)=centers(m,1)-round(LX(m,1)*10);
            centers(m,2)=centers(m,2)-round(LY(m,1)*10)-7;
            continue;
        end
    end
    
    centers(centers(:,1)<=65,:)=[];
    centers(centers(:,1)>=295,:)=[];
    
    flushdata(vid);
    
    d=[];
    dT=[];
    ED=cell(1,1);
    JPnum=cell(1,1);
    Jnum=[];
    
    for m=1:size(centers,1)
        for n=1:size(centers,1) 
             d(n,1)=sqrt((centers(n,1)-centers(m,1))^2+(centers(n,2)-centers(m,2))^2);
        end
        ED{m,1}=d;   
    end
    
    for m=1:size(centers,1)
        JPnum{m,1}=find(ED{m,1}<20&ED{m,1}>0);
    end

    for m=1:size(centers,1)
        Jnum(m,1)=size(JPnum{m,1},1);
    end
    netnum=zeros(size(centers,1),1);
    
    for i=1:size(centers,1)
        if netnum(i,1)==0
            netnum(i,1)=max(netnum)+1;
        else
            netnum(i,1)=netnum(i,1);
        end
        for j=1:size(JPnum{i,1},1)
            netnum(JPnum{i,1}(j,1),1)=netnum(i,1);
        end
    end
    
    Target=[0,360];
    TargetChange=[360,360];
    
    for m=1:size(centers,1)   
        dT(m,1)=abs(sqrt((Target(1,1)-centers(m,1))^2+(Target(1,2)-centers(m,2))^2));
    end
    
    netCA=[netnum,Jnum,centers,dT];
    netCB=sortrows(netCA,[1 5]);
    netCC=[];
    for mn=2:size(netCB)
        netCC(1,1)=1;
        if netCB(mn,1)==netCB(mn-1,1)
            netCC(mn,1)=netCC(mn-1,1);    
        else
            netCC(mn,1)=netCC(mn-1,1)+1;
        end
    end
    netCB(:,1)=netCC;
    
    for in=1:size(netCB(:,1))
        dTC=[];
        netcenter=sum(netCB(netCB(:,1)==in,3))/sum(netCB(:,1)==in);
        netCBX=netCB(netCB(:,1)==in,3);
        netCBY=netCB(netCB(:,1)==in,4);
        if netcenter>=120 && netcenter<=163 || netcenter>=196 && netcenter<=239
            for ic=1:sum(netCB(:,1)==in)
                dTC(ic,1)=abs(sqrt((TargetChange(1,1)-netCBX(ic))^2+(TargetChange(1,2)-netCBY(ic))^2)); 
            end
            netCB(netCB(:,1)==in,5)=dTC;
        else
            continue
        end
    end
    
    netCB=sortrows(netCB,[1 5]);
    
    conum=[];
    cop=[];
    copl=[];
    for mm=1:max(netCB(:,1))
        co=[];
        co=find(netCB(:,1)==mm);  
        [di,Li]=max(netCB(co,5));
        conum(mm,1)=size(co,1);
        for nn=round(conum(mm,1)/2)+1:conum(mm,1)
            copt=[];
            if mm==1
                if nn==0
                    continue
                else
                copt(nn,1)=nn;
                end
            else 
                copt(nn,1)=nn+sum(conum(1:mm-1,1));
            end        
            cop=[cop;copt];   
         end
    end
    cop(cop==0)=[];
    
    for mn=1:size(cop,1)
        copl(mn,:)=netCB(cop(mn,1),3:4);
    end
    
    R=16;
    for x=1:R
        for y=1:R
            D = sqrt( abs(x-round(R/2))^2 + abs(y-round(R/2))^2 );
            if( D<=R/2  )
                   C(x,y)=255;
            else
                C(x,y,:)=0;
            end
        end
    end
    
    centerMO=copl;
    for i=1:size(centerMO,1)
        for j=centerMO(i,2)-ceil(R/2)+1:centerMO(i,2)+floor(R/2)
            for k=centerMO(i,1)-ceil(R/2)+1:centerMO(i,1)+floor(R/2)
                if k<=0
                    k=1;
                else
                    k=k;
                end
                P(j,k,2:3)=P(j,k,2:3)+0;
                P(j,k,:)=P(j,k,1)+C(j-(centerMO(i,2)-ceil(R/2)+1)+1,k-(centerMO(i,1)-ceil(R/2)+1)+1);
            end
        end
    end
    
    set(MAINFIG,'CData',P);
    refreshdata(MAINFIG);
    drawnow limitrate 
    pause(0.4);
    
    P=STphoto*gl2;
    set(MAINFIG,'CData',P);
    refreshdata(MAINFIG);
    pause(0.3);
       
    toc
end

drawnow