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

STphoto=imread('Interactive light field basemap.png');
axis off
P=ones(360,360,3);
gl2=0.2;
P=STphoto(5:355,5:255,:)*gl2;

MAINFIG=image(P);
axis off
drawnow;
pause
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
    gl1=0;
    LSD=0.5;
    P=STphoto*gl2;
    
    P_ini=getsnapshot(vid);
    P_cut=P_ini(15:1065,460:1510,:);
    P_rgb=rgb2gray(P_cut);
    P_TW = im2bw(P_rgb,0.7);
    P_cut(:,:,1)=0;
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
    
    flushdata(vid);
    
    Target=sum(centers)/size(centers,1);
    
    d=[];
    dT=[];
    ED=cell(1,1);
    JPnum=cell(1,1);
    Jnum=[];
    for m=1:size(centers,1)
        dT(m,1)=sqrt((Target(1,1)-centers(m,1))^2+(Target(1,2)-centers(m,2))^2);
    end
    
    for m=1:size(centers,1)
        for n=1:size(centers,1) 
             d(n,1)=sqrt((centers(n,1)-centers(m,1))^2+(centers(n,2)-centers(m,2))^2);
        end
        ED{m,1}=d;   
    end
    
    for m=1:size(centers,1)
        JPnum{m,1}=find(ED{m,1}<18&ED{m,1}>0);
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
    
    Pt=ones(360,360,3)*0;
    Dr=9;
    for x=1:Dr
        for y=1:Dr
            D = sqrt( abs(x-round(Dr/2))^2 + abs(y-round(Dr/2))^2 );
            if( D<=Dr/2  )
                C(x,y)=255;
            else
                C(x,y,:)=0;
            end
        end
    end
    
    centerMt=netCB(:,3:4);
    for in=size(centerMt(:,1))
        if centerMt(in,1)>=360-Dr-2 || centerMt(in,1)<=Dr || centerMt(in,2)>=360-Dr-2 || centerMt(in,2)<=Dr
            centerMt(in,1:2)=0;       
        else
            centerMt(in,1:2)=centerMt(in,1:2);
        end
    end
    centerMt(centerMt(:,1)==0,:)=[];
    
    for i=1:size(centerMt,1)
        for j=centerMt(i,2)-ceil(Dr/2)+1:centerMt(i,2)+floor(Dr/2)
            for k=centerMt(i,1)-ceil(Dr/2)+1:centerMt(i,1)+floor(Dr/2)
                Pt(j,k,2:3)=Pt(j,k,2:3)+0;
                Pt(j,k,:)=Pt(j,k,1)+C(j-(centerMt(i,2)-ceil(Dr/2)+1)+1,k-(centerMt(i,1)-ceil(Dr/2)+1)+1);
            end
        end
    end
    
    R=14;
    
    figure(10);
    imshow(Pt);
    [xm,ym]=ginput(1);
    for im=size(xm(:,1))
        if xm(im,1)>=360-R
           xm(im,1)=360-R;
        elseif xm(im,1)<=R
            xm(im,1)=0+R;
        else
            xm(im,1)=xm(im,1);
        end

        if ym(im,1)>=360-R
           ym(im,1)=360-R;
        elseif ym(im,1)<=R
           ym(im,1)=0+R;
        else
           ym(im,1)=ym(im,1);
        end
    end
        
    XYm=[xm,ym];
    
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
    
    centerMO=round(XYm);
    for i=1:size(centerMO,1)
        for j=centerMO(i,2)-ceil(R/2)+1:centerMO(i,2)+floor(R/2)
            for k=centerMO(i,1)-ceil(R/2)+1:centerMO(i,1)+floor(R/2)
                P(j,k,2:3)=P(j,k,2:3)+0;
                P(j,k,:)=P(j,k,1)+C(j-(centerMO(i,2)-ceil(R/2)+1)+1,k-(centerMO(i,1)-ceil(R/2)+1)+1);
            end
        end
    end
    
    set(MAINFIG,'CData',P);
    refreshdata(MAINFIG);
    drawnow limitrate 
    pause(0.3);
    
    P=STphoto*gl2;
    set(MAINFIG,'CData',P);
    drawnow
    refreshdata(MAINFIG);
    
    toc
end

drawnow