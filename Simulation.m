%Kinetics Simulation of Magbot
%Changeable Uneven Light Field & Nonequilibrium Drive
%Extended Attract Region
clear all
close all
clc
tic
rng('shuffle');

%Parameter Input
DataSavePath='D:\SimulationData\Magbot\';
FileName='Magbot_KUDE100_420_23071409';
ShowMode='off';

StepTimes=0.004;
SubStepTimes=0.004;
SubStepNumber=StepTimes/SubStepTimes;
RecordTimeIntervals=0.2;
RecordTimeIntervalStep=RecordTimeIntervals/StepTimes;
RecordId=0;
TotalTimes=15*60;
TotalTimeStep=TotalTimes/StepTimes+1;

SpatialScalemm=1.25;
Lengthmm=420;
LengthPixel=Lengthmm/SpatialScalemm;
RadiusPixel=LengthPixel/2;
PositionValue=(1:LengthPixel)-0.5;
[GridCenterX,GridCenterY]=meshgrid(PositionValue,PositionValue);
Distance=sqrt((GridCenterX-RadiusPixel).^2+(GridCenterY-RadiusPixel).^2);
XDistance=LengthPixel-GridCenterX;
YDistance=LengthPixel-GridCenterY;
Environment=double(Distance<=RadiusPixel);
%Environment=ones(LengthPixel,LengthPixel);
Scene=255*Environment;
BoundaryMode=1;

MagbotNumber=100;
MagbotRadiusmm=10;
MagbotRadiusPixel=MagbotRadiusmm/SpatialScalemm;
PolarityType=3;
PolarityAngle=2*pi/PolarityType;
RotateMatrix=[cos(PolarityAngle),-sin(PolarityAngle);sin(PolarityAngle),cos(PolarityAngle)];
ActingLimit=0;

RotationBasalSpeedDistribution=ones(LengthPixel,LengthPixel);%Distance/RadiusPixel;%XDistance/LengthPixel;%
RotationBasalSpeedrpsMax=5;
RotationBasalSpeedrps=RotationBasalSpeedrpsMax*RotationBasalSpeedDistribution.*Environment; %Corresponding to Light Intensity
RotationBasalSpeedradpstep=RotationBasalSpeedrps*2*pi*StepTimes;
RotationSpeedCal=@(mOmega,R,N)0.5*mOmega*((R/MagbotRadiusPixel).^(-1.2)+N.^(-0.7));
RotationNoiseMaxRatio=0.05;

TranslationalNoiseMaxSpeedRatio=1;
TranslationalSigma=25/StepTimes;

AttractRatio=2; %Corresponding to magnetic strength
AttractDistanceCmm=AttractRatio*1.5;
StrongCouplingDistanceCmm=AttractRatio*0.6;
StrongCouplingAnglerad=asin(0.05*StrongCouplingDistanceCmm/MagbotRadiusmm)*2;

BreakLevel=10^(-1.2*AttractRatio);
DifferenceC=10;

InvadeLimit=0.5;

SwitchTimes=30*60;
SwitchTimeStep=SwitchTimes/StepTimes;

ErrorToleration=1e-4;

MagbotRecord=cell(1,TotalTimes/RecordTimeIntervals+1);
BondRecord=cell(1,TotalTimes/RecordTimeIntervals+1);
PolymerRecord=cell(1,TotalTimes/RecordTimeIntervals+1);
InputParameter=cell(2,26);
InputParameter(1,:)={'StepTimes','SubStepTimes','RecordTimeIntervals','TotalTimes','SpatialScalemm','Lengthmm','BoundaryMode','MagbotNumber','MagbotRadiusmm','PolarityType','ActingLimit','RotationBasalSpeedDistribution','RotationBasalSpeedrpsMax','RotationBasalSpeedrps','RotationSpeedCal','RotationNoiseMaxRatio','TranslationalNoiseMaxSpeedRatio','TranslationalSigma','AttractRatio','AttractDistanceCmm','StrongCouplingDistanceCmm','BreakLevel','DifferenceC','InvadeLimit','SwitchTimes','ErrorToleration'};
InputParameter(2,:)={StepTimes,SubStepTimes,RecordTimeIntervals,TotalTimes,SpatialScalemm,Lengthmm,BoundaryMode,MagbotNumber,MagbotRadiusmm,PolarityType,ActingLimit,RotationBasalSpeedDistribution,RotationBasalSpeedrpsMax,RotationBasalSpeedrps,RotationSpeedCal,RotationNoiseMaxRatio,TranslationalNoiseMaxSpeedRatio,TranslationalSigma,AttractRatio,AttractDistanceCmm,StrongCouplingDistanceCmm,BreakLevel,DifferenceC,InvadeLimit,SwitchTimes,ErrorToleration};

CirclePointNumber=20;
CirclePointAngle=2*pi*linspace(0,1,CirclePointNumber+1);
ScalebarPixel=round(100/SpatialScalemm);
ScaleBarPointX=[LengthPixel-7-ScalebarPixel,LengthPixel-7,LengthPixel-7,LengthPixel-7-ScalebarPixel];
ScaleBarPointY=[LengthPixel-20,LengthPixel-20,LengthPixel-25,LengthPixel-25];
% MovieRobots=VideoWriter(fullfile(DataSavePath,[FileName,'Movie']),'Motion JPEG AVI');
% MovieRobots.Quality=75;
% MovieRobots.FrameRate = 20;

%Magbot Initialization
PolymerMass=ones(1,MagbotNumber);
load('MatrixInitialPosition100_420.mat')
MagbotStatus=cell(MagbotNumber,1);
for n=1:MagbotNumber
    MagbotStatus{n,1}{1,1}=zeros(1+PolarityType,11); %[id,absolute x,absolute y,relative x,relative y,radius,energy,connected bond number,nearest neighbor number,rspeed,actingC;connection state,absolute x,absolute y,relative x,relative y,nearest neighbor id,nearest bond id,gap distance,bond energy,rspeed difference,average rspeed]
    MagbotStatus{n,1}{2,1}=1; %individual number
    MagbotStatus{n,1}{2,2}=InitialPosition(n,:)/SpatialScalemm; %rotation center
    MagbotStatus{n,1}{2,3}=MagbotRadiusPixel; %rotation radius
    MagbotStatus{n,1}{2,4}=RotationBasalSpeedradpstep(ceil(MagbotStatus{n,1}{2,2}(2)),ceil(MagbotStatus{n,1}{2,2}(1))); %final rotation speed
    MagbotStatus{n,1}{2,5}=0.25*MagbotStatus{n,1}{2,3}^2*MagbotStatus{n,1}{2,4}^2; %Total Energy
    MagbotStatus{n,1}{2,6}=[rand*0.5*pi+(randi(40)-21)*0.5*pi,0]; %Translational Angle, Persistent Step
    MagbotStatus{n,1}{2,7}=zeros(MagbotStatus{n,1}{2,1},MagbotStatus{n,1}{2,1}); %Connection Matrix
    MagbotStatus{n,1}{2,8}=0; %Connected Bond Number
    MagbotStatus{n,1}{2,9}=n; %Mapping Relationship
    MagbotStatus{n,1}{2,10}=MagbotStatus{n,1}{2,4}; %Average rotation speed
    MagbotStatus{n,1}{2,11}=[0,0]; %Directional Move
    
    MagbotStatus{n,1}{1,1}(1,1)=n;
    MagbotStatus{n,1}{1,1}(1,2:3)=InitialPosition(n,:)/SpatialScalemm;
    MagbotStatus{n,1}{1,1}(1,4:6)=[0,0,0];
    MagbotStatus{n,1}{1,1}(1,7)=0.25*MagbotStatus{n,1}{2,3}^2*MagbotStatus{n,1}{2,4}^2;
    MagbotStatus{n,1}{1,1}(1,8)=0;
    MagbotStatus{n,1}{1,1}(1,9)=0;
    MagbotStatus{n,1}{1,1}(1,10)=MagbotStatus{n,1}{2,4};
    MagbotStatus{n,1}{1,1}(1,11)=ActingLimit+(1-ActingLimit)*rand;
    
    BondAngle=rand*2*pi;
    for p=1:PolarityType
        MagbotStatus{n,1}{1,1}(1+p,1)=0;
        MagbotStatus{n,1}{1,1}(1+p,2:3)=MagbotStatus{n,1}{1,1}(1,2:3)+MagbotRadiusPixel*[cos(BondAngle),sin(BondAngle)];
        MagbotStatus{n,1}{1,1}(1+p,4:5)=MagbotRadiusPixel*[cos(BondAngle),sin(BondAngle)];
        MagbotStatus{n,1}{1,1}(1+p,6:11)=0;
        BondAngle=BondAngle+PolarityAngle;
    end
end

%Run
% open(MovieRobots);
% figure('visible',ShowMode);
for i_t=1:TotalTimeStep
    
    %Light Field Change
    if rem(i_t,SwitchTimeStep)==0
        RotationBasalSpeedDistribution=1-RotationBasalSpeedDistribution;
        RotationBasalSpeedrps=RotationBasalSpeedrpsMax*RotationBasalSpeedDistribution.*Environment;
        RotationBasalSpeedradpstep=RotationBasalSpeedrps*2*pi*StepTimes;
        for i_p=1:size(MagbotStatus,1)
            rspeedsum=0;
            dmove=[0,0];
            dmovewight=0;
            for i_n=1:MagbotStatus{i_p,1}{2,1}
                MagbotStatus{i_p,1}{1,i_n}(1,10)=RotationBasalSpeedradpstep(ceil(MagbotStatus{i_p,1}{1,i_n}(1,3)),ceil(MagbotStatus{i_p,1}{1,i_n}(1,2)));
                rspeedsum=rspeedsum+MagbotStatus{i_p,1}{1,i_n}(1,10);
                dmove=dmove+MagbotStatus{i_p,1}{1,i_n}(1,11)*MagbotStatus{i_p,1}{1,i_n}(1,6)*MagbotStatus{i_p,1}{1,i_n}(1,10)*(MagbotStatus{i_p,1}{2,2}-[MagbotStatus{i_p,1}{1,i_n}(1,2),MagbotStatus{i_p,1}{1,i_n}(1,3)]);
                dmovewight=dmovewight+MagbotStatus{i_p,1}{1,i_n}(1,11)*MagbotStatus{i_p,1}{1,i_n}(1,6)^2*MagbotStatus{i_p,1}{1,i_n}(1,10);
            end
            MagbotStatus{i_p,1}{2,10}=rspeedsum/MagbotStatus{i_p,1}{2,1};
            MagbotStatus{i_p,1}{2,11}=[0,0];
            if dmovewight>ErrorToleration&&norm(dmove)>ErrorToleration
                dmove=dmove/dmovewight;
                MagbotStatus{i_p,1}{2,11}=dmove/norm(dmove);
            end
            MagbotStatus{i_p,1}{2,4}=RotationSpeedCal(MagbotStatus{i_p,1}{2,10},MagbotStatus{i_p,1}{2,3},MagbotStatus{i_p,1}{2,1});
            MagbotStatus{i_p,1}{2,5}=0;
            for i_n=1:MagbotStatus{i_p,1}{2,1}
                R=max(MagbotStatus{i_p,1}{1,i_n}(1,6),sqrt(2)/2);
                MagbotStatus{i_p,1}{1,i_n}(1,7)=0.5*(R*MagbotStatus{i_p,1}{2,4})^2;
                MagbotStatus{i_p,1}{2,5}=MagbotStatus{i_p,1}{2,5}+MagbotStatus{i_p,1}{1,i_n}(1,7);
            end
            for i_n=1:MagbotStatus{i_p,1}{2,1}
                for i_l=1:PolarityType
                    if MagbotStatus{i_p,1}{1,i_n}(1+i_l,1)==1
                        other=find(MagbotStatus{i_p,1}{2,9}==MagbotStatus{i_p,1}{1,i_n}(1+i_l,6));
                        MagbotStatus{i_p,1}{1,i_n}(1+i_l,9)=(MagbotStatus{i_p,1}{1,i_n}(1,7)+MagbotStatus{i_p,1}{1,other}(1,7))/(MagbotStatus{i_p,1}{1,i_n}(1,8)+MagbotStatus{i_p,1}{1,other}(1,8)-1);
                        MagbotStatus{i_p,1}{1,i_n}(1+i_l,10)=abs(MagbotStatus{i_p,1}{1,i_n}(1,10)-MagbotStatus{i_p,1}{1,other}(1,10));
                        MagbotStatus{i_p,1}{1,i_n}(1+i_l,11)=(MagbotStatus{i_p,1}{1,i_n}(1,10)+MagbotStatus{i_p,1}{1,other}(1,10))/2;
                    end
                end
            end
        end
    end
    
    %Record
    if rem(i_t-1,RecordTimeIntervalStep)==0
        MagbotAll=zeros(MagbotNumber,10);
        BondAll=zeros(PolarityType*MagbotNumber,11);
        PolymerConection=zeros(MagbotNumber,MagbotNumber);
        for i_p=1:size(MagbotStatus,1)
            for i_n=1:MagbotStatus{i_p,1}{2,1}
                MagbotAll(MagbotStatus{i_p,1}{1,i_n}(1,1),:)=MagbotStatus{i_p,1}{1,i_n}(1,2:11);
                for i_m=1:PolarityType
                    BondAll(i_m+PolarityType*(MagbotStatus{i_p,1}{1,i_n}(1,1)-1),:)=MagbotStatus{i_p,1}{1,i_n}(1+i_m,:);
                    LinkState=MagbotStatus{i_p,1}{1,i_n}(1+i_m,1);
                    if LinkState==1
                        PolymerConection(MagbotStatus{i_p,1}{1,i_n}(1,1),MagbotStatus{i_p,1}{1,i_n}(1+i_m,6))=1;
                    end
                end
            end
        end
        RecordId=RecordId+1;
        MagbotRecord{1,RecordId}=MagbotAll;
        BondRecord{1,RecordId}=BondAll;
        PolymerRecord{1,RecordId}=PolymerConection;
        
        %Show
%         TimeElapsed=(i_t-1)*StepTimes/60/60/24;
%         imshow(uint8(Scene.*RotationBasalSpeedDistribution));
%         daspect([1 1 1]);
%         hold on
%         for i_n=1:MagbotNumber
%             RobotEdgeX=MagbotAll(i_n,1)+MagbotRadiusPixel*cos(CirclePointAngle);
%             RobotEdgeY=MagbotAll(i_n,2)+MagbotRadiusPixel*sin(CirclePointAngle);
%             fill(RobotEdgeX,RobotEdgeY,[0 0 1]);
%         end
%         for i_b=1:PolarityType*MagbotNumber
%             plot(BondAll(i_b,2),BondAll(i_b,3),'.','MarkerEdgeColor',(1-BondAll(i_b,1))*[0 1 0]+BondAll(i_b,1)*[1 0 0],'MarkerSize',5);
%         end
%         TimeStr=strcat('Real Time:',32,datestr(TimeElapsed,'HH:MM:SS.FFF'));
%         text(2,12,TimeStr,'Color',[1,0,0],'FontSize',15);
%         fill(ScaleBarPointX,ScaleBarPointY,[1,0,0]);
%         ScaleStr='100 mm';
%         text(LengthPixel-ScalebarPixel-32,LengthPixel-10,ScaleStr,'Color',[1,0,0],'FontSize',15);
%         drawnow;
%         Frame=getframe(gcf);
%         writeVideo(MovieRobots,Frame);
%         clf

i_t
    end
    
    %Rotation and Random Move
    RotationMove=zeros(size(MagbotStatus,1),1);
    TranslationalMove=zeros(size(MagbotStatus,1),2);
    for i_n=1:size(MagbotStatus,1)
        RotationMove(i_n)=(1+randsrc*rand*RotationNoiseMaxRatio)*MagbotStatus{i_n,1}{2,4};
        if MagbotStatus{i_n,1}{2,6}(2)<=0
            MagbotStatus{i_n,1}{2,6}(1)=rand*0.5*pi+(randi(40)-21)*0.5*pi;
            MagbotStatus{i_n,1}{2,6}(2)=ceil(abs(normrnd(0,TranslationalSigma)))-1;
        else
            MagbotStatus{i_n,1}{2,6}(2)=MagbotStatus{i_n,1}{2,6}(2)-1;
        end
        TranslationalLength=TranslationalNoiseMaxSpeedRatio*MagbotStatus{i_n,1}{2,10}/(2*pi)/SpatialScalemm;
        RandomLength=1/MagbotStatus{i_n,1}{2,1}*TranslationalLength*rand;
        DirectionalLength=TranslationalLength-RandomLength;
        TranslationalMove(i_n,1)=RandomLength*cos(MagbotStatus{i_n,1}{2,6}(1))+DirectionalLength*MagbotStatus{i_n,1}{2,11}(1);
        TranslationalMove(i_n,2)=RandomLength*sin(MagbotStatus{i_n,1}{2,6}(1))+DirectionalLength*MagbotStatus{i_n,1}{2,11}(2);
    end
    
    for i_ts=1:SubStepNumber
        for i_n=1:size(MagbotStatus,1)
            MagbotStatus{i_n,1}{2,2}=MagbotStatus{i_n,1}{2,2}+TranslationalMove(i_n,:)/SubStepNumber;
            [NewCenterAll,~]=meshgrid(MagbotStatus{i_n,1}{2,2},1:(1+PolarityType));
            RotateMatrixTemp=[cos(RotationMove(i_n)/SubStepNumber),-sin(RotationMove(i_n)/SubStepNumber);sin(RotationMove(i_n)/SubStepNumber),cos(RotationMove(i_n)/SubStepNumber)];
            OldCenter=zeros(MagbotStatus{i_n,1}{2,1},2);
            for i_m=1:MagbotStatus{i_n,1}{2,1}
                for i_l=1:1+PolarityType
                    MagbotStatus{i_n,1}{1,i_m}(i_l,4:5)=(RotateMatrixTemp*MagbotStatus{i_n,1}{1,i_m}(i_l,4:5)')';
                end
                MagbotStatus{i_n,1}{1,i_m}(:,2:3)=MagbotStatus{i_n,1}{1,i_m}(:,4:5)+NewCenterAll;
                OldCenter(i_m,:)=MagbotStatus{i_n,1}{1,i_m}(1,2:3);
            end
            
            %Collide
            if MagbotStatus{i_n,1}{2,1}>1
                for i_nn=2:MagbotStatus{i_n,1}{2,1}
                    for i_nm=(i_nn-1):-1:1
                        distance=norm(MagbotStatus{i_n,1}{1,i_nn}(1,2:3)-MagbotStatus{i_n,1}{1,i_nm}(1,2:3));
                        if distance<2*MagbotRadiusPixel-ErrorToleration
                            BackDistance=2*MagbotRadiusPixel-distance;
                            if distance>ErrorToleration
                                MagbotStatus{i_n,1}{1,i_nn}(1,2:3)=MagbotStatus{i_n,1}{1,i_nn}(1,2:3)+(MagbotStatus{i_n,1}{1,i_nn}(1,2:3)-MagbotStatus{i_n,1}{1,i_nm}(1,2:3))/distance*BackDistance;
                            else
                                RandomDirection=rand*2*pi;
                                MagbotStatus{i_n,1}{1,i_nn}(1,2:3)=MagbotStatus{i_n,1}{1,i_nn}(1,2:3)+BackDistance*[cos(RandomDirection),sin(RandomDirection)];
                            end
                        end
                    end
                end
            end
            
            if i_n>1
                for i_nn=1:MagbotStatus{i_n,1}{2,1}
                    for i_m=i_n-1:-1:1
                        for i_mn=1:MagbotStatus{i_m,1}{2,1}
                            distance=norm(MagbotStatus{i_n,1}{1,i_nn}(1,2:3)-MagbotStatus{i_m,1}{1,i_mn}(1,2:3));
                            if distance<2*MagbotRadiusPixel-ErrorToleration
                                BackDistance=2*MagbotRadiusPixel-distance;
                                if distance>ErrorToleration
                                    BackVector=(MagbotStatus{i_n,1}{1,i_nn}(1,2:3)-MagbotStatus{i_m,1}{1,i_mn}(1,2:3))/distance*BackDistance;
                                    for i_k=1:MagbotStatus{i_n,1}{2,1}
                                        MagbotStatus{i_n,1}{1,i_k}(1,2:3)=MagbotStatus{i_n,1}{1,i_k}(1,2:3)+BackVector;
                                    end
                                else
                                    RandomDirection=rand*2*pi;
                                    for i_k=1:MagbotStatus{i_n,1}{2,1}
                                        MagbotStatus{i_n,1}{1,i_k}(1,2:3)=MagbotStatus{i_n,1}{1,i_k}(1,2:3)+BackDistance*[cos(RandomDirection),sin(RandomDirection)];
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            for i_nn=1:MagbotStatus{i_n,1}{2,1}
                if BoundaryMode==1
                    CenterDistance=norm(MagbotStatus{i_n,1}{1,i_nn}(1,2:3)-[RadiusPixel,RadiusPixel]);
                    if (CenterDistance+MagbotRadiusPixel)>RadiusPixel+ErrorToleration
                        BackDistance=(CenterDistance+MagbotRadiusPixel)-RadiusPixel;
                        BackVector=([RadiusPixel,RadiusPixel]-MagbotStatus{i_n,1}{1,i_nn}(1,2:3))/CenterDistance*BackDistance;
                        for i_k=1:MagbotStatus{i_n,1}{2,1}
                            MagbotStatus{i_n,1}{1,i_k}(1,2:3)=MagbotStatus{i_n,1}{1,i_k}(1,2:3)+BackVector;
                        end
                    end
                elseif BoundaryMode==2
                    LeftDistance=MagbotStatus{i_n,1}{1,i_nn}(1,2);
                    if LeftDistance<MagbotRadiusPixel-ErrorToleration
                        BackDistance=MagbotRadiusPixel-LeftDistance;
                        for i_k=1:MagbotStatus{i_n,1}{2,1}
                            MagbotStatus{i_n,1}{1,i_k}(1,2)=MagbotStatus{i_n,1}{1,i_k}(1,2)+BackDistance;
                        end
                    end
                    RightDistance=LengthPixel-MagbotStatus{i_n,1}{1,i_nn}(1,2);
                    if RightDistance<MagbotRadiusPixel-ErrorToleration
                        BackDistance=MagbotRadiusPixel-RightDistance;
                        for i_k=1:MagbotStatus{i_n,1}{2,1}
                            MagbotStatus{i_n,1}{1,i_k}(1,2)=MagbotStatus{i_n,1}{1,i_k}(1,2)-BackDistance;
                        end
                    end
                    UpDistance=MagbotStatus{i_n,1}{1,i_nn}(1,3);
                    if UpDistance<MagbotRadiusPixel-ErrorToleration
                        BackDistance=MagbotRadiusPixel-UpDistance;
                        for i_k=1:MagbotStatus{i_n,1}{2,1}
                            MagbotStatus{i_n,1}{1,i_k}(1,3)=MagbotStatus{i_n,1}{1,i_k}(1,3)+BackDistance;
                        end
                    end
                    DownDistance=LengthPixel-MagbotStatus{i_n,1}{1,i_nn}(1,3);
                    if DownDistance<MagbotRadiusPixel-ErrorToleration
                        BackDistance=MagbotRadiusPixel-DownDistance;
                        for i_k=1:MagbotStatus{i_n,1}{2,1}
                            MagbotStatus{i_n,1}{1,i_k}(1,3)=MagbotStatus{i_n,1}{1,i_k}(1,3)-BackDistance;
                        end
                    end
                end
            end
            
            MagbotStatus{i_n,1}{2,2}=0;
            for i_nn=1:MagbotStatus{i_n,1}{2,1}
                MagbotStatus{i_n,1}{2,2}=MagbotStatus{i_n,1}{2,2}+MagbotStatus{i_n,1}{1,i_nn}(1,2:3);
                for i_l=1:PolarityType
                    MagbotStatus{i_n,1}{1,i_nn}(i_l+1,2:3)=MagbotStatus{i_n,1}{1,i_nn}(i_l+1,2:3)+MagbotStatus{i_n,1}{1,i_nn}(1,2:3)-OldCenter(i_nn,:);
                end
            end
            MagbotStatus{i_n,1}{2,2}=MagbotStatus{i_n,1}{2,2}/MagbotStatus{i_n,1}{2,1};
            for i_m=1:MagbotStatus{i_n,1}{2,1}
                for i_l=1:1+PolarityType
                    MagbotStatus{i_n,1}{1,i_m}(i_l,4:5)=MagbotStatus{i_n,1}{1,i_m}(i_l,2:3)-MagbotStatus{i_n,1}{2,2};
                end
            end
        end
    end
    
    %Weak Attract
    for i_n=1:size(MagbotStatus,1)
        if MagbotStatus{i_n,1}{2,1}>1
            for i_nn=1:MagbotStatus{i_n,1}{2,1}
                AttractNumber=0;
                AttractPoint=[0,0];
                AttractBodyPoint=[0,0];
                AttractGapDistance=LengthPixel*ones(PolarityType,2);
                for i_nm=1:MagbotStatus{i_n,1}{2,1}
                    if i_nm~=i_nn
                        distance=norm(MagbotStatus{i_n,1}{1,i_nn}(1,2:3)-MagbotStatus{i_n,1}{1,i_nm}(1,2:3));
                        if distance<(2*MagbotRadiusPixel+AttractDistanceCmm/SpatialScalemm-ErrorToleration)
                            Bonddistance=zeros(PolarityType,PolarityType);
                            for i_k=1:PolarityType
                                for i_l=1:PolarityType
                                    Bonddistance(i_k,i_l)=norm(MagbotStatus{i_n,1}{1,i_nn}(1+i_k,2:3)-MagbotStatus{i_n,1}{1,i_nm}(1+i_l,2:3));
                                end
                            end
                            GapDistance=min(Bonddistance(:));
                            [MinIdn,MinIdm]=find(Bonddistance==GapDistance,1);
                            if GapDistance<AttractDistanceCmm/SpatialScalemm-ErrorToleration
                                AttractNumber=AttractNumber+1;
                                AttractPoint(AttractNumber,1:2)=MagbotStatus{i_n,1}{1,i_nm}(1+MinIdm,2:3);
                                AttractBodyPoint(AttractNumber,1:2)=MagbotStatus{i_n,1}{1,i_nm}(1,2:3);
                                AttractGapDistance(AttractNumber,:)=[GapDistance,MinIdn];
                            end
                        end
                    end
                end
                if AttractNumber>0
                    [MinGap,MinGapId]=min(AttractGapDistance(:,1));
                    [OwnAngle,~]=cart2pol(MagbotStatus{i_n,1}{1,i_nn}(1+AttractGapDistance(MinGapId,2),2)-MagbotStatus{i_n,1}{1,i_nn}(1,2),MagbotStatus{i_n,1}{1,i_nn}(1+AttractGapDistance(MinGapId,2),3)-MagbotStatus{i_n,1}{1,i_nn}(1,3));
                    OwnPointX=MagbotStatus{i_n,1}{1,i_nn}(1,2)+MagbotRadiusPixel*cos([OwnAngle-StrongCouplingAnglerad,OwnAngle,OwnAngle+StrongCouplingAnglerad]);
                    OwnPointY=MagbotStatus{i_n,1}{1,i_nn}(1,3)+MagbotRadiusPixel*sin([OwnAngle-StrongCouplingAnglerad,OwnAngle,OwnAngle+StrongCouplingAnglerad]);
                    [OtherAngle,~]=cart2pol(AttractPoint(MinGapId,1)-AttractBodyPoint(MinGapId,1),AttractPoint(MinGapId,2)-AttractBodyPoint(MinGapId,2));
                    OtherPointX=AttractBodyPoint(MinGapId,1)+MagbotRadiusPixel*cos([OtherAngle-StrongCouplingAnglerad,OtherAngle,OtherAngle+StrongCouplingAnglerad]);
                    OtherPointY=AttractBodyPoint(MinGapId,2)+MagbotRadiusPixel*sin([OtherAngle-StrongCouplingAnglerad,OtherAngle,OtherAngle+StrongCouplingAnglerad]);
                    DGapdistance=zeros(3,3);
                    for i_k=1:3
                        for i_l=1:3
                            DGapdistance(i_k,i_l)=norm([OwnPointX(i_k)-OtherPointX(i_l),OwnPointY(i_k)-OtherPointY(i_l)]);
                        end
                    end
                    MinDGapdistance=min(DGapdistance(:));
                    [MinIdn,MinIdm]=find(DGapdistance==MinDGapdistance,1);
                    MagbotStatus{i_n,1}{1,i_nn}(1,2:3)=2*[OtherPointX(MinIdm),OtherPointY(MinIdm)]-AttractBodyPoint(MinGapId,1:2);
                    ContactDirection=([OtherPointX(MinIdm),OtherPointY(MinIdm)]-MagbotStatus{i_n,1}{1,i_nn}(1,2:3))/norm([OtherPointX(MinIdm),OtherPointY(MinIdm)]-MagbotStatus{i_n,1}{1,i_nn}(1,2:3));
                    [ContactAngle,~]=cart2pol(ContactDirection(1),ContactDirection(2));
                    if MinIdn==1
                        BondAngle=ContactAngle+StrongCouplingAnglerad;
                    elseif MinIdn==2
                        BondAngle=ContactAngle;
                    elseif MinIdn==3
                        BondAngle=ContactAngle-StrongCouplingAnglerad;
                    end
                    MagbotStatus{i_n,1}{1,i_nn}(2,2:3)=MagbotStatus{i_n,1}{1,i_nn}(1,2:3)+MagbotRadiusPixel*[cos(BondAngle),sin(BondAngle)];
                    for i_l=2:PolarityType
                        BondAngle=BondAngle+PolarityAngle;
                        MagbotStatus{i_n,1}{1,i_nn}(1+i_l,2:3)=MagbotStatus{i_n,1}{1,i_nn}(1,2:3)+MagbotRadiusPixel*[cos(BondAngle),sin(BondAngle)];
                    end
                end
            end
        end
    end
    
    if size(MagbotStatus,1)>1
        for i_n=2:size(MagbotStatus,1)
            if MagbotStatus{i_n,1}{2,1}>1
                for i_nn=1:MagbotStatus{i_n,1}{2,1}
                    for i_m=i_n-1:-1:1
                        for i_mn=1:MagbotStatus{i_m,1}{2,1}
                            distance=norm(MagbotStatus{i_n,1}{1,i_nn}(1,2:3)-MagbotStatus{i_m,1}{1,i_mn}(1,2:3));
                            if distance<(2*MagbotRadiusPixel+AttractDistanceCmm/SpatialScalemm-ErrorToleration)
                                Bonddistance=zeros(PolarityType,PolarityType);
                                for i_k=1:PolarityType
                                    for i_l=1:PolarityType
                                        Bonddistance(i_k,i_l)=norm(MagbotStatus{i_n,1}{1,i_nn}(1+i_k,2:3)-MagbotStatus{i_m,1}{1,i_mn}(1+i_l,2:3));
                                    end
                                end
                                GapDistance=min(Bonddistance(:));
                                [MinIdn,MinIdm]=find(Bonddistance==GapDistance,1);
                                if GapDistance<AttractDistanceCmm/SpatialScalemm-ErrorToleration
                                    AttractVector=2*MagbotStatus{i_m,1}{1,i_mn}(1+MinIdm,2:3)-MagbotStatus{i_m,1}{1,i_mn}(1,2:3)-MagbotStatus{i_n,1}{1,i_nn}(1,2:3);
                                    for i_k=1:MagbotStatus{i_n,1}{2,1}
                                        for i_l=1:1+PolarityType
                                            MagbotStatus{i_n,1}{1,i_k}(i_l,2:3)=MagbotStatus{i_n,1}{1,i_k}(i_l,2:3)+AttractVector/MagbotStatus{i_n,1}{2,1};
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            else
                i_nn=1;
                AttractNumber=0;
                AttractPoint=[0,0];
                AttractBodyPoint=[0,0];
                AttractGapDistance=LengthPixel*ones(PolarityType,2);
                for i_m=i_n-1:-1:1
                    for i_mn=1:MagbotStatus{i_m,1}{2,1}
                        distance=norm(MagbotStatus{i_n,1}{1,i_nn}(1,2:3)-MagbotStatus{i_m,1}{1,i_mn}(1,2:3));
                        if distance<(2*MagbotRadiusPixel+AttractDistanceCmm/SpatialScalemm-ErrorToleration)
                            Bonddistance=zeros(PolarityType,PolarityType);
                            for i_k=1:PolarityType
                                for i_l=1:PolarityType
                                    Bonddistance(i_k,i_l)=norm(MagbotStatus{i_n,1}{1,i_nn}(1+i_k,2:3)-MagbotStatus{i_m,1}{1,i_mn}(1+i_l,2:3));
                                end
                            end
                            GapDistance=min(Bonddistance(:));
                            [MinIdn,MinIdm]=find(Bonddistance==GapDistance,1);
                            if GapDistance<AttractDistanceCmm/SpatialScalemm-ErrorToleration
                                AttractNumber=AttractNumber+1;
                                AttractPoint(AttractNumber,1:2)=MagbotStatus{i_m,1}{1,i_mn}(1+MinIdm,2:3);
                                AttractBodyPoint(AttractNumber,1:2)=MagbotStatus{i_m,1}{1,i_mn}(1,2:3);
                                AttractGapDistance(AttractNumber,:)=[GapDistance,MinIdn];
                            end
                        end
                    end
                end
                if AttractNumber>0
                    [MinGap,MinGapId]=min(AttractGapDistance(:,1));
                    [OwnAngle,~]=cart2pol(MagbotStatus{i_n,1}{1,i_nn}(1+AttractGapDistance(MinGapId,2),2)-MagbotStatus{i_n,1}{1,i_nn}(1,2),MagbotStatus{i_n,1}{1,i_nn}(1+AttractGapDistance(MinGapId,2),3)-MagbotStatus{i_n,1}{1,i_nn}(1,3));
                    OwnPointX=MagbotStatus{i_n,1}{1,i_nn}(1,2)+MagbotRadiusPixel*cos([OwnAngle-StrongCouplingAnglerad,OwnAngle,OwnAngle+StrongCouplingAnglerad]);
                    OwnPointY=MagbotStatus{i_n,1}{1,i_nn}(1,3)+MagbotRadiusPixel*sin([OwnAngle-StrongCouplingAnglerad,OwnAngle,OwnAngle+StrongCouplingAnglerad]);
                    [OtherAngle,~]=cart2pol(AttractPoint(MinGapId,1)-AttractBodyPoint(MinGapId,1),AttractPoint(MinGapId,2)-AttractBodyPoint(MinGapId,2));
                    OtherPointX=AttractBodyPoint(MinGapId,1)+MagbotRadiusPixel*cos([OtherAngle-StrongCouplingAnglerad,OtherAngle,OtherAngle+StrongCouplingAnglerad]);
                    OtherPointY=AttractBodyPoint(MinGapId,2)+MagbotRadiusPixel*sin([OtherAngle-StrongCouplingAnglerad,OtherAngle,OtherAngle+StrongCouplingAnglerad]);
                    DGapdistance=zeros(3,3);
                    for i_k=1:3
                        for i_l=1:3
                            DGapdistance(i_k,i_l)=norm([OwnPointX(i_k)-OtherPointX(i_l),OwnPointY(i_k)-OtherPointY(i_l)]);
                        end
                    end
                    MinDGapdistance=min(DGapdistance(:));
                    [MinIdn,MinIdm]=find(DGapdistance==MinDGapdistance,1);
                    MagbotStatus{i_n,1}{1,i_nn}(1,2:3)=2*[OtherPointX(MinIdm),OtherPointY(MinIdm)]-AttractBodyPoint(MinGapId,1:2);
                    ContactDirection=([OtherPointX(MinIdm),OtherPointY(MinIdm)]-MagbotStatus{i_n,1}{1,i_nn}(1,2:3))/norm([OtherPointX(MinIdm),OtherPointY(MinIdm)]-MagbotStatus{i_n,1}{1,i_nn}(1,2:3));
                    [ContactAngle,~]=cart2pol(ContactDirection(1),ContactDirection(2));
                    if MinIdn==1
                        BondAngle=ContactAngle+StrongCouplingAnglerad;
                    elseif MinIdn==2
                        BondAngle=ContactAngle;
                    elseif MinIdn==3
                        BondAngle=ContactAngle-StrongCouplingAnglerad;
                    end
                    MagbotStatus{i_n,1}{1,i_nn}(2,2:3)=MagbotStatus{i_n,1}{1,i_nn}(1,2:3)+MagbotRadiusPixel*[cos(BondAngle),sin(BondAngle)];
                    for i_l=2:PolarityType
                        BondAngle=BondAngle+PolarityAngle;
                        MagbotStatus{i_n,1}{1,i_nn}(1+i_l,2:3)=MagbotStatus{i_n,1}{1,i_nn}(1,2:3)+MagbotRadiusPixel*[cos(BondAngle),sin(BondAngle)];
                    end
                end
            end
        end
    end
    
    %Collide
    for i_n=1:size(MagbotStatus,1)
        OldCenter=zeros(MagbotStatus{i_n,1}{2,1},2);
        for i_m=1:MagbotStatus{i_n,1}{2,1}
            OldCenter(i_m,:)=MagbotStatus{i_n,1}{1,i_m}(1,2:3);
        end
        
        if MagbotStatus{i_n,1}{2,1}>1
            for i_nn=2:MagbotStatus{i_n,1}{2,1}
                for i_nm=(i_nn-1):-1:1
                    distance=norm(MagbotStatus{i_n,1}{1,i_nn}(1,2:3)-MagbotStatus{i_n,1}{1,i_nm}(1,2:3));
                    if distance<2*MagbotRadiusPixel-ErrorToleration
                        BackDistance=2*MagbotRadiusPixel-distance;
                        if distance>ErrorToleration
                            MagbotStatus{i_n,1}{1,i_nn}(1,2:3)=MagbotStatus{i_n,1}{1,i_nn}(1,2:3)+(MagbotStatus{i_n,1}{1,i_nn}(1,2:3)-MagbotStatus{i_n,1}{1,i_nm}(1,2:3))/distance*BackDistance;
                        else
                            RandomDirection=rand*2*pi;
                            MagbotStatus{i_n,1}{1,i_nn}(1,2:3)=MagbotStatus{i_n,1}{1,i_nn}(1,2:3)+BackDistance*[cos(RandomDirection),sin(RandomDirection)];
                        end
                    end
                end
            end
        end
        
        if i_n>1
            for i_nn=1:MagbotStatus{i_n,1}{2,1}
                for i_m=i_n-1:-1:1
                    for i_mn=1:MagbotStatus{i_m,1}{2,1}
                        distance=norm(MagbotStatus{i_n,1}{1,i_nn}(1,2:3)-MagbotStatus{i_m,1}{1,i_mn}(1,2:3));
                        if distance<2*MagbotRadiusPixel-ErrorToleration
                            BackDistance=2*MagbotRadiusPixel-distance;
                            if distance>ErrorToleration
                                BackVector=(MagbotStatus{i_n,1}{1,i_nn}(1,2:3)-MagbotStatus{i_m,1}{1,i_mn}(1,2:3))/distance*BackDistance;
                                for i_k=1:MagbotStatus{i_n,1}{2,1}
                                    MagbotStatus{i_n,1}{1,i_k}(1,2:3)=MagbotStatus{i_n,1}{1,i_k}(1,2:3)+BackVector;
                                end
                            else
                                RandomDirection=rand*2*pi;
                                for i_k=1:MagbotStatus{i_n,1}{2,1}
                                    MagbotStatus{i_n,1}{1,i_k}(1,2:3)=MagbotStatus{i_n,1}{1,i_k}(1,2:3)+BackDistance*[cos(RandomDirection),sin(RandomDirection)];
                                end
                            end
                        end
                    end
                end
            end
        end
        
        for i_nn=1:MagbotStatus{i_n,1}{2,1}
            if BoundaryMode==1
                CenterDistance=norm(MagbotStatus{i_n,1}{1,i_nn}(1,2:3)-[RadiusPixel,RadiusPixel]);
                if (CenterDistance+MagbotRadiusPixel)>RadiusPixel+ErrorToleration
                    BackDistance=(CenterDistance+MagbotRadiusPixel)-RadiusPixel;
                    BackVector=([RadiusPixel,RadiusPixel]-MagbotStatus{i_n,1}{1,i_nn}(1,2:3))/CenterDistance*BackDistance;
                    for i_k=1:MagbotStatus{i_n,1}{2,1}
                        MagbotStatus{i_n,1}{1,i_k}(1,2:3)=MagbotStatus{i_n,1}{1,i_k}(1,2:3)+BackVector;
                    end
                end
            elseif BoundaryMode==2
                LeftDistance=MagbotStatus{i_n,1}{1,i_nn}(1,2);
                if LeftDistance<MagbotRadiusPixel-ErrorToleration
                    BackDistance=MagbotRadiusPixel-LeftDistance;
                    for i_k=1:MagbotStatus{i_n,1}{2,1}
                        MagbotStatus{i_n,1}{1,i_k}(1,2)=MagbotStatus{i_n,1}{1,i_k}(1,2)+BackDistance;
                    end
                end
                RightDistance=LengthPixel-MagbotStatus{i_n,1}{1,i_nn}(1,2);
                if RightDistance<MagbotRadiusPixel-ErrorToleration
                    BackDistance=MagbotRadiusPixel-RightDistance;
                    for i_k=1:MagbotStatus{i_n,1}{2,1}
                        MagbotStatus{i_n,1}{1,i_k}(1,2)=MagbotStatus{i_n,1}{1,i_k}(1,2)-BackDistance;
                    end
                end
                UpDistance=MagbotStatus{i_n,1}{1,i_nn}(1,3);
                if UpDistance<MagbotRadiusPixel-ErrorToleration
                    BackDistance=MagbotRadiusPixel-UpDistance;
                    for i_k=1:MagbotStatus{i_n,1}{2,1}
                        MagbotStatus{i_n,1}{1,i_k}(1,3)=MagbotStatus{i_n,1}{1,i_k}(1,3)+BackDistance;
                    end
                end
                DownDistance=LengthPixel-MagbotStatus{i_n,1}{1,i_nn}(1,3);
                if DownDistance<MagbotRadiusPixel-ErrorToleration
                    BackDistance=MagbotRadiusPixel-DownDistance;
                    for i_k=1:MagbotStatus{i_n,1}{2,1}
                        MagbotStatus{i_n,1}{1,i_k}(1,3)=MagbotStatus{i_n,1}{1,i_k}(1,3)-BackDistance;
                    end
                end
            end
        end
        
        MagbotStatus{i_n,1}{2,2}=0;
        for i_nn=1:MagbotStatus{i_n,1}{2,1}
            MagbotStatus{i_n,1}{2,2}=MagbotStatus{i_n,1}{2,2}+MagbotStatus{i_n,1}{1,i_nn}(1,2:3);
            for i_l=1:PolarityType
                MagbotStatus{i_n,1}{1,i_nn}(i_l+1,2:3)=MagbotStatus{i_n,1}{1,i_nn}(i_l+1,2:3)-OldCenter(i_nn,:)+MagbotStatus{i_n,1}{1,i_nn}(1,2:3);
            end
        end
        MagbotStatus{i_n,1}{2,2}=MagbotStatus{i_n,1}{2,2}/MagbotStatus{i_n,1}{2,1};
        for i_m=1:MagbotStatus{i_n,1}{2,1}
            for i_l=1:1+PolarityType
                MagbotStatus{i_n,1}{1,i_m}(i_l,4:5)=MagbotStatus{i_n,1}{1,i_m}(i_l,2:3)-MagbotStatus{i_n,1}{2,2};
            end
        end
    end
    
    NewMonomer=zeros(MagbotNumber,14);
    NewBondPosition=zeros(MagbotNumber*PolarityType,2);
    for i_n=1:size(MagbotStatus,1)
        for i_m=1:MagbotStatus{i_n,1}{2,1}
            NewMonomer(MagbotStatus{i_n,1}{1,i_m}(1,1),1)=i_n;
            NewMonomer(MagbotStatus{i_n,1}{1,i_m}(1,1),2:11)=MagbotStatus{i_n,1}{1,i_m}(1,2:11);
            NewMonomer(MagbotStatus{i_n,1}{1,i_m}(1,1),12)=MagbotStatus{i_n,1}{2,1};
            NewMonomer(MagbotStatus{i_n,1}{1,i_m}(1,1),13)=MagbotStatus{i_n,1}{2,6}(1)+MagbotStatus{i_n,1}{2,1}*(randi(3)-2)*2*pi;
            NewMonomer(MagbotStatus{i_n,1}{1,i_m}(1,1),14)=MagbotStatus{i_n,1}{2,6}(2);
            for i_l=1:PolarityType
                NewBondPosition(i_l+PolarityType*(MagbotStatus{i_n,1}{1,i_m}(1,1)-1),:)=MagbotStatus{i_n,1}{1,i_m}(1+i_l,2:3);
            end
        end
    end
    
    InvadeAll=zeros(MagbotNumber,MagbotNumber);
    for i_n=1:MagbotNumber
        for i_m=i_n:MagbotNumber
            if i_m~=i_n
                InvadeAll(i_n,i_m)=(2*MagbotRadiusPixel-norm(NewMonomer(i_n,2:3)-NewMonomer(i_m,2:3)))/MagbotRadiusPixel;
                InvadeAll(i_m,i_n)=InvadeAll(i_n,i_m);
            else
                if BoundaryMode==1
                    InvadeAll(i_n,i_m)=(norm(NewMonomer(i_n,2:3)-[RadiusPixel,RadiusPixel])+MagbotRadiusPixel-RadiusPixel)/MagbotRadiusPixel;
                elseif BoundaryMode==2
                    LeftDistance=MagbotRadiusPixel-NewMonomer(i_n,2);
                    RightDistance=MagbotRadiusPixel-(LengthPixel-NewMonomer(i_n,2));
                    UpDistance=MagbotRadiusPixel-NewMonomer(i_n,3);
                    DownDistance=MagbotRadiusPixel-(LengthPixel-NewMonomer(i_n,3));
                    InvadeAll(i_n,i_m)=max([LeftDistance,RightDistance,UpDistance,DownDistance])/MagbotRadiusPixel;
                end
            end
        end
    end
    
    %Deterministic bonding according to distance
    ConnectionMatrixBig=zeros(MagbotNumber,MagbotNumber);
    NearestN=zeros(MagbotNumber,1);
    Gap=LengthPixel*ones(PolarityType*MagbotNumber,1);
    for i_n=1:(MagbotNumber-1)
        for i_m=(i_n+1):MagbotNumber
            distance=norm(NewMonomer(i_n,2:3)-NewMonomer(i_m,2:3));
            if (distance-2*MagbotRadiusPixel)<StrongCouplingDistanceCmm/SpatialScalemm-ErrorToleration
                NearestN(i_n)=NearestN(i_n)+1;
                NearestN(i_m)=NearestN(i_m)+1;
                if max(InvadeAll(i_n,:))<InvadeLimit-ErrorToleration&&max(InvadeAll(i_m,:))<InvadeLimit-ErrorToleration
                    Bonddistance=zeros(PolarityType,PolarityType);
                    for i_k=1:PolarityType
                        for i_l=1:PolarityType
                            Bonddistance(i_k,i_l)=norm(NewBondPosition(i_k+PolarityType*(i_n-1),:)-NewBondPosition(i_l+PolarityType*(i_m-1),:));
                        end
                    end
                    GapDistance=min(Bonddistance(:));
                    [MinIdn,MinIdm]=find(Bonddistance==GapDistance,1);
                    if GapDistance<StrongCouplingDistanceCmm/SpatialScalemm-ErrorToleration
                        ConnectionMatrixBig(i_n,i_m)=MinIdn;
                        ConnectionMatrixBig(i_m,i_n)=MinIdm;
                        Gap(MinIdn+PolarityType*(i_n-1))=GapDistance;
                        Gap(MinIdm+PolarityType*(i_m-1))=GapDistance;
                    end
                end
            end
        end
    end
    
    ConnectionMatrixBigS=double(ConnectionMatrixBig>0);
    ConnectionGraph=graph(ConnectionMatrixBigS);
    ConnectionState=conncomp(ConnectionGraph);
    PolymerNumber=max(ConnectionState);
    MagbotStatus=cell(PolymerNumber,1);
    for i_n=1:PolymerNumber
        MagbotStatus{i_n,1}{2,1}=sum(double(ConnectionState==i_n));
        MappingRelationship=zeros(1,MagbotStatus{i_n,1}{2,1});
        Inid=0;
        PositionSum=[0,0];
        rspeedsum=0;
        dmove=[0,0];
        dmovewight=0;
        MagbotStatus{i_n,1}{2,6}=[0,0];
        for i_m=1:MagbotNumber
            if ConnectionState(i_m)==i_n
                Inid=Inid+1;
                MappingRelationship(Inid)=i_m;
                MagbotStatus{i_n,1}{1,Inid}=zeros(1+PolarityType,11);
                MagbotStatus{i_n,1}{1,Inid}(1,1)=i_m;
                MagbotStatus{i_n,1}{1,Inid}(1,2:3)=NewMonomer(i_m,2:3);
                PositionSum=PositionSum+NewMonomer(i_m,2:3);
                MagbotStatus{i_n,1}{2,6}=MagbotStatus{i_n,1}{2,6}+NewMonomer(i_m,13:14);
                MagbotStatus{i_n,1}{1,Inid}(1,8)=sum(ConnectionMatrixBigS(i_m,:));
                MagbotStatus{i_n,1}{1,Inid}(1,9)=NearestN(i_m);
                MagbotStatus{i_n,1}{1,Inid}(1,10)=RotationBasalSpeedradpstep(ceil(MagbotStatus{i_n,1}{1,Inid}(1,3)),ceil(MagbotStatus{i_n,1}{1,Inid}(1,2)));
                MagbotStatus{i_n,1}{1,Inid}(1,11)=NewMonomer(i_m,11);
                rspeedsum=rspeedsum+MagbotStatus{i_n,1}{1,Inid}(1,10);
                for i_l=1:PolarityType
                    MagbotStatus{i_n,1}{1,Inid}(1+i_l,2:3)=NewBondPosition(i_l+PolarityType*(i_m-1),:);
                    MagbotStatus{i_n,1}{1,Inid}(1+i_l,8)=Gap(i_l+PolarityType*(i_m-1));
                end
                for i_k=1:MagbotNumber
                    if ConnectionMatrixBig(i_m,i_k)>0
                        MagbotStatus{i_n,1}{1,Inid}(1+ConnectionMatrixBig(i_m,i_k),1)=1;
                        MagbotStatus{i_n,1}{1,Inid}(1+ConnectionMatrixBig(i_m,i_k),6)=i_k;
                        MagbotStatus{i_n,1}{1,Inid}(1+ConnectionMatrixBig(i_m,i_k),7)=ConnectionMatrixBig(i_k,i_m);
                    end
                end
            end
        end
        MagbotStatus{i_n,1}{2,10}=rspeedsum/MagbotStatus{i_n,1}{2,1};
        MagbotStatus{i_n,1}{2,2}=PositionSum/MagbotStatus{i_n,1}{2,1};
        MagbotStatus{i_n,1}{2,6}=MagbotStatus{i_n,1}{2,6}/MagbotStatus{i_n,1}{2,1};
        Radius=zeros(1,MagbotStatus{i_n,1}{2,1});
        for i_m=1:MagbotStatus{i_n,1}{2,1}
            Radius(i_m)=norm(MagbotStatus{i_n,1}{1,i_m}(1,2:3)-MagbotStatus{i_n,1}{2,2});
            for i_l=1:1+PolarityType
                MagbotStatus{i_n,1}{1,i_m}(i_l,4:5)=MagbotStatus{i_n,1}{1,i_m}(i_l,2:3)-MagbotStatus{i_n,1}{2,2};
            end
            MagbotStatus{i_n,1}{1,i_m}(1,6)=Radius(i_m);
            dmove=dmove+MagbotStatus{i_n,1}{1,i_m}(1,11)*MagbotStatus{i_n,1}{1,i_m}(1,6)*MagbotStatus{i_n,1}{1,i_m}(1,10)*(MagbotStatus{i_n,1}{2,2}-[MagbotStatus{i_n,1}{1,i_m}(1,2),MagbotStatus{i_n,1}{1,i_m}(1,3)]);
            dmovewight=dmovewight+MagbotStatus{i_n,1}{1,i_m}(1,11)*MagbotStatus{i_n,1}{1,i_m}(1,6)^2*MagbotStatus{i_n,1}{1,i_m}(1,10);
        end
        [NewR,NewO]=sort(Radius,'ascend');
        NewMagbotStatus=cell(1,MagbotStatus{i_n,1}{2,1});
        for i_nn=1:MagbotStatus{i_n,1}{2,1}
            NewMagbotStatus{1,i_nn}=MagbotStatus{i_n,1}{1,NewO(i_nn)};
        end
        for i_nn=1:MagbotStatus{i_n,1}{2,1}
            MagbotStatus{i_n,1}{1,i_nn}=NewMagbotStatus{1,i_nn};
        end
        NewMappingRelationship=MappingRelationship(NewO);
        MagbotStatus{i_n,1}{2,9}=NewMappingRelationship;
        MagbotStatus{i_n,1}{2,11}=[0,0];
        if dmovewight>ErrorToleration&&norm(dmove)>ErrorToleration
            dmove=dmove/dmovewight;
            MagbotStatus{i_n,1}{2,11}=dmove/norm(dmove);
        end
        MagbotStatus{i_n,1}{2,3}=max(Radius)+MagbotRadiusPixel;
        MagbotStatus{i_n,1}{2,4}=RotationSpeedCal(MagbotStatus{i_n,1}{2,10},MagbotStatus{i_n,1}{2,3},MagbotStatus{i_n,1}{2,1});
        MagbotStatus{i_n,1}{2,5}=0;
        for i_m=1:MagbotStatus{i_n,1}{2,1}
            R=max(MagbotStatus{i_n,1}{1,i_m}(1,6),sqrt(2)/2);
            MagbotStatus{i_n,1}{1,i_m}(1,7)=0.5*(R*MagbotStatus{i_n,1}{2,4})^2;
            MagbotStatus{i_n,1}{2,5}=MagbotStatus{i_n,1}{2,5}+MagbotStatus{i_n,1}{1,i_m}(1,7);
        end
        MagbotStatus{i_n,1}{2,7}=zeros(MagbotStatus{i_n,1}{2,1},MagbotStatus{i_n,1}{2,1});
        for i_m=1:MagbotStatus{i_n,1}{2,1}
            for i_l=1:PolarityType
                if MagbotStatus{i_n,1}{1,i_m}(1+i_l,1)==1
                    other=find(NewMappingRelationship==MagbotStatus{i_n,1}{1,i_m}(1+i_l,6));
                    MagbotStatus{i_n,1}{1,i_m}(1+i_l,9)=(MagbotStatus{i_n,1}{1,i_m}(1,7)+MagbotStatus{i_n,1}{1,other}(1,7))/(MagbotStatus{i_n,1}{1,i_m}(1,8)+MagbotStatus{i_n,1}{1,other}(1,8)-1);
                    MagbotStatus{i_n,1}{1,i_m}(1+i_l,10)=abs(MagbotStatus{i_n,1}{1,i_m}(1,10)-MagbotStatus{i_n,1}{1,other}(1,10));
                    MagbotStatus{i_n,1}{1,i_m}(1+i_l,11)=(MagbotStatus{i_n,1}{1,i_m}(1,10)+MagbotStatus{i_n,1}{1,other}(1,10))/2;
                end
            end
            for i_k=1:MagbotStatus{i_n,1}{2,1}
                MagbotStatus{i_n,1}{2,7}(i_m,i_k)=ConnectionMatrixBigS(NewMappingRelationship(i_m),NewMappingRelationship(i_k));
            end
        end
        Link=MagbotStatus{i_n,1}{2,7};
        MagbotStatus{i_n,1}{2,8}=sum(Link(:))/2;
    end
    
    %Probabilistic breaking according to energy
    for i_p=size(MagbotStatus,1):-1:1
        for i_n=1:MagbotStatus{i_p,1}{2,1}
            for i_b=1:PolarityType
                if MagbotStatus{i_p,1}{1,i_n}(1+i_b,1)>0
                    if rand<((BreakLevel*(MagbotStatus{i_p,1}{1,i_n}(1+i_b,11)+DifferenceC*MagbotStatus{i_p,1}{1,i_n}(1+i_b,10)))-ErrorToleration)
                        MagbotStatus{i_p,1}{1,i_n}(1+i_b,1)=0;
                        MagbotStatus{i_p,1}{1,i_n}(1+i_b,9)=0;
                        MagbotStatus{i_p,1}{1,i_n}(1+i_b,10)=0;
                        MagbotStatus{i_p,1}{1,i_n}(1+i_b,11)=0;
                        MagbotStatus{i_p,1}{1,i_n}(1,8)=MagbotStatus{i_p,1}{1,i_n}(1,8)-1;
                        other=find(MagbotStatus{i_p,1}{2,9}==MagbotStatus{i_p,1}{1,i_n}(1+i_b,6));
                        MagbotStatus{i_p,1}{1,other}(1+MagbotStatus{i_p,1}{1,i_n}(1+i_b,7),1)=0;
                        MagbotStatus{i_p,1}{1,other}(1+MagbotStatus{i_p,1}{1,i_n}(1+i_b,7),9)=0;
                        MagbotStatus{i_p,1}{1,other}(1+MagbotStatus{i_p,1}{1,i_n}(1+i_b,7),10)=0;
                        MagbotStatus{i_p,1}{1,other}(1+MagbotStatus{i_p,1}{1,i_n}(1+i_b,7),11)=0;
                        MagbotStatus{i_p,1}{1,other}(1,8)=MagbotStatus{i_p,1}{1,other}(1,8)-1;
                        MagbotStatus{i_p,1}{2,7}(i_n,other)=0;
                        MagbotStatus{i_p,1}{2,7}(other,i_n)=0;
                        MagbotStatus{i_p,1}{2,8}=MagbotStatus{i_p,1}{2,8}-1;
                    end
                end
            end
        end
        ConnectionGraph=graph(MagbotStatus{i_p,1}{2,7});
        ConnectionState=conncomp(ConnectionGraph);
        PolymerNumber=max(ConnectionState);
        if PolymerNumber>1
            NewPolymer=cell(PolymerNumber,1);
            for i_pp=1:PolymerNumber
                NewPolymer{i_pp,1}{2,1}=sum(double(ConnectionState==i_pp));
                MappingRelationship=zeros(1,NewPolymer{i_pp,1}{2,1});
                MappingRelationshipS=zeros(1,NewPolymer{i_pp,1}{2,1});
                Inid=0;
                PositionSum=[0,0];
                rspeedsum=0;
                dmove=[0,0];
                dmovewight=0;
                for i_m=1:MagbotStatus{i_p,1}{2,1}
                    if ConnectionState(i_m)==i_pp
                        Inid=Inid+1;
                        MappingRelationship(Inid)=MagbotStatus{i_p,1}{2,9}(1,i_m);
                        MappingRelationshipS(Inid)=i_m;
                        NewPolymer{i_pp,1}{1,Inid}=zeros(1+PolarityType,11);
                        NewPolymer{i_pp,1}{1,Inid}(1,1)=MagbotStatus{i_p,1}{2,9}(1,i_m);
                        NewPolymer{i_pp,1}{1,Inid}(1,2:3)=NewMonomer(MagbotStatus{i_p,1}{2,9}(1,i_m),2:3);
                        PositionSum=PositionSum+NewMonomer(MagbotStatus{i_p,1}{2,9}(1,i_m),2:3);
                        NewPolymer{i_pp,1}{1,Inid}(1,8)=sum(MagbotStatus{i_p,1}{2,7}(i_m,:));
                        NewPolymer{i_pp,1}{1,Inid}(1,9)=NearestN(MagbotStatus{i_p,1}{2,9}(1,i_m));
                        NewPolymer{i_pp,1}{1,Inid}(1,10)=RotationBasalSpeedradpstep(ceil(NewPolymer{i_pp,1}{1,Inid}(1,3)),ceil(NewPolymer{i_pp,1}{1,Inid}(1,2)));
                        NewPolymer{i_pp,1}{1,Inid}(1,11)=NewMonomer(MagbotStatus{i_p,1}{2,9}(1,i_m),11);
                        rspeedsum=rspeedsum+NewPolymer{i_pp,1}{1,Inid}(1,10);
                        for i_l=1:PolarityType
                            NewPolymer{i_pp,1}{1,Inid}(1+i_l,2:3)=NewBondPosition(i_l+PolarityType*(MagbotStatus{i_p,1}{2,9}(1,i_m)-1),:);
                            NewPolymer{i_pp,1}{1,Inid}(1+i_l,8)=Gap(i_l+PolarityType*(MagbotStatus{i_p,1}{2,9}(1,i_m)-1));
                        end
                        for i_k=1:MagbotStatus{i_p,1}{2,1}
                            if MagbotStatus{i_p,1}{2,7}(i_m,i_k)>0
                                NewPolymer{i_pp,1}{1,Inid}(1+ConnectionMatrixBig(MagbotStatus{i_p,1}{2,9}(1,i_m),MagbotStatus{i_p,1}{2,9}(1,i_k)),1)=1;
                                NewPolymer{i_pp,1}{1,Inid}(1+ConnectionMatrixBig(MagbotStatus{i_p,1}{2,9}(1,i_m),MagbotStatus{i_p,1}{2,9}(1,i_k)),6)=MagbotStatus{i_p,1}{2,9}(1,i_k);
                                NewPolymer{i_pp,1}{1,Inid}(1+ConnectionMatrixBig(MagbotStatus{i_p,1}{2,9}(1,i_m),MagbotStatus{i_p,1}{2,9}(1,i_k)),7)=ConnectionMatrixBig(MagbotStatus{i_p,1}{2,9}(1,i_k),MagbotStatus{i_p,1}{2,9}(1,i_m));
                            end
                        end
                    end
                end
                NewPolymer{i_pp,1}{2,10}=rspeedsum/NewPolymer{i_pp,1}{2,1};
                NewPolymer{i_pp,1}{2,2}=PositionSum/NewPolymer{i_pp,1}{2,1};
                NewPolymer{i_pp,1}{2,6}=[MagbotStatus{i_p,1}{2,6}(1)+(randi(3)-2)*2*pi,ceil(abs(normrnd(0,TranslationalSigma)))];
                Radius=zeros(1,NewPolymer{i_pp,1}{2,1});
                for i_m=1:NewPolymer{i_pp,1}{2,1}
                    Radius(i_m)=norm(NewPolymer{i_pp,1}{1,i_m}(1,2:3)-NewPolymer{i_pp,1}{2,2});
                    for i_l=1:1+PolarityType
                        NewPolymer{i_pp,1}{1,i_m}(i_l,4:5)=NewPolymer{i_pp,1}{1,i_m}(i_l,2:3)-NewPolymer{i_pp,1}{2,2};
                    end
                    NewPolymer{i_pp,1}{1,i_m}(1,6)=Radius(i_m);
                    dmove=dmove+NewPolymer{i_pp,1}{1,i_m}(1,11)*NewPolymer{i_pp,1}{1,i_m}(1,6)*NewPolymer{i_pp,1}{1,i_m}(1,10)*(NewPolymer{i_pp,1}{2,2}-[NewPolymer{i_pp,1}{1,i_m}(1,2),NewPolymer{i_pp,1}{1,i_m}(1,3)]);
                    dmovewight=dmovewight+NewPolymer{i_pp,1}{1,i_m}(1,11)*NewPolymer{i_pp,1}{1,i_m}(1,6)^2*NewPolymer{i_pp,1}{1,i_m}(1,10);
                end
                [NewR,NewO]=sort(Radius,'ascend');
                NewMagbotStatus=cell(1,NewPolymer{i_pp,1}{2,1});
                for i_nn=1:NewPolymer{i_pp,1}{2,1}
                    NewMagbotStatus{1,i_nn}=NewPolymer{i_pp,1}{1,NewO(i_nn)};
                end
                for i_nn=1:NewPolymer{i_pp,1}{2,1}
                    NewPolymer{i_pp,1}{1,i_nn}=NewMagbotStatus{1,i_nn};
                end
                NewMappingRelationship=MappingRelationship(NewO);
                NewMappingRelationshipS=MappingRelationshipS(NewO);
                NewPolymer{i_pp,1}{2,9}=NewMappingRelationship;
                NewPolymer{i_pp,1}{2,11}=[0,0];
                if dmovewight>ErrorToleration&&norm(dmove)>ErrorToleration
                    dmove=dmove/dmovewight;
                    NewPolymer{i_pp,1}{2,11}=dmove/norm(dmove);
                end
                NewPolymer{i_pp,1}{2,3}=max(Radius)+MagbotRadiusPixel;
                NewPolymer{i_pp,1}{2,4}=RotationSpeedCal(NewPolymer{i_pp,1}{2,10},NewPolymer{i_pp,1}{2,3},NewPolymer{i_pp,1}{2,1});
                NewPolymer{i_pp,1}{2,5}=0;
                for i_m=1:NewPolymer{i_pp,1}{2,1}
                    R=max(NewPolymer{i_pp,1}{1,i_m}(1,6),sqrt(2)/2);
                    NewPolymer{i_pp,1}{1,i_m}(1,7)=0.5*(R*NewPolymer{i_pp,1}{2,4})^2;
                    NewPolymer{i_pp,1}{2,5}=NewPolymer{i_pp,1}{2,5}+NewPolymer{i_pp,1}{1,i_m}(1,7);
                end
                NewPolymer{i_pp,1}{2,7}=zeros(NewPolymer{i_pp,1}{2,1},NewPolymer{i_pp,1}{2,1});
                for i_m=1:NewPolymer{i_pp,1}{2,1}
                    for i_l=1:PolarityType
                        if NewPolymer{i_pp,1}{1,i_m}(1+i_l,1)==1
                            other=find(NewMappingRelationship==NewPolymer{i_pp,1}{1,i_m}(1+i_l,6));
                            NewPolymer{i_pp,1}{1,i_m}(1+i_l,9)=(NewPolymer{i_pp,1}{1,i_m}(1,7)+NewPolymer{i_pp,1}{1,other}(1,7))/(NewPolymer{i_pp,1}{1,i_m}(1,8)+NewPolymer{i_pp,1}{1,other}(1,8)-1);
                            NewPolymer{i_pp,1}{1,i_m}(1+i_l,10)=abs(NewPolymer{i_pp,1}{1,i_m}(1,10)-NewPolymer{i_pp,1}{1,other}(1,10));
                            NewPolymer{i_pp,1}{1,i_m}(1+i_l,11)=(NewPolymer{i_pp,1}{1,i_m}(1,10)+NewPolymer{i_pp,1}{1,other}(1,10))/2;
                        end
                    end
                    for i_k=1:NewPolymer{i_pp,1}{2,1}
                        NewPolymer{i_pp,1}{2,7}(i_m,i_k)=MagbotStatus{i_p,1}{2,7}(NewMappingRelationshipS(i_m),NewMappingRelationshipS(i_k));
                    end
                end
                Link=NewPolymer{i_pp,1}{2,7};
                NewPolymer{i_pp,1}{2,8}=sum(Link(:))/2;
                MagbotStatus{size(MagbotStatus,1)+1,1}=NewPolymer{i_pp,1};
            end
            MagbotStatus(i_p)=[];
        end
    end
    
    PolymerMass=ones(1,size(MagbotStatus,1));
    for i_n=1:size(MagbotStatus,1)
        PolymerMass(i_n)=MagbotStatus{i_n,1}{2,1};
    end
    [NewPolymerMass,NewOrder]=sort(PolymerMass,'descend');
    NewMagbotStatus=cell(size(MagbotStatus,1),1);
    for i_n=1:size(MagbotStatus,1)
        NewMagbotStatus{i_n,1}=MagbotStatus{NewOrder(i_n),1};
    end
    MagbotStatus=NewMagbotStatus;
    
    
end

% close(MovieRobots)
% figure('visible','on');
% close all

Results=cell(2,3);
Results(1,:)={'MagbotRecord','BondRecord','PolymerRecord'};
Results(2,:)={MagbotRecord,BondRecord,PolymerRecord};

ExportData=cell(2,1);
ExportData{1}=InputParameter;
ExportData{2}=Results;

filename=[DataSavePath,FileName,'Data.mat'];
save(filename,'ExportData');
toc