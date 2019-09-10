%real data plot
clear all
close all

yy=importdata('C:\Users\Arka\Desktop\pso_ulti_final\real_data\Residual gravity anomaly of the offshore Louisiana salt dome_observed data.bln');
x_obs1=yy(:,1);
data_grv=yy(:,2);
x_obs=linspace(-6110.11236731062,6836.48208131634,length(yy));
data1=spline(x_obs1',data_grv',x_obs);
z_obs=0;
data1=data1*10^-5;
data2=diff(data1(:))./diff(x_obs(:));
data2=data2*10^8;
data3=spline(x_obs(1:end-1),data2,x_obs);
data_gv=data1*10^5;
data_gr=data3;
figure(1)
plot(x_obs1,data_grv)
hold on
plot(x_obs,data_gv)
figure(2)
plot(x_obs,data_gr)
data_all=[data_gv, data_gr];
for plg=4:7
for cnt=1:65
    cnt
            %% Problem Definition
CostFunction =@(x,data) myCostFunction(x,data)+...
    10000*(constrained1(x)+constrained2(x)+constrained3(x)+constrained4(x)+constrained5(x)+constrained6(x)+...
          constrained7(x)+constrained8(x));

%nVar=8;                %Number of Unknown Variable 
nVar=plg*2+1;                 %Number of Unknown Variable including density

VarSize = [1 nVar];     %Matrix size of Decision variables

VarMin= -ones(1,nVar);          %Lower Bound of Unknown variable
VarMax= ones(1,nVar);           %Upper Bound of Unknown variable

%% Parameters of PSO

MaxIt = 300;           %Maximum number of iterations

nPoP =  50;            %Population size or swarm size


w=1;                    %inertia coefficient

c1=1;                   %Personal Accelaration 
c2=2;                   %Social Accelaration

%% Initialization

Empty.Particle.Position =[];
Empty.Particle.Velocity =[];
Empty.Particle.Cost     =[];

Empty.Particle.Best.Position =[];
Empty.Particle.Best.Cost     =[];

Particle=repmat(Empty.Particle, nPoP,1);

%initial global best
GlobalBest.Cost= Inf;

for i=1:nPoP
    %initialize position with random number from VarMin and VarMax
    for j=1:nVar
        Particle(i).Position(j) =(VarMax(j)-VarMin(j))*rand(1) + VarMin(j);
    end
    
    %Initialize Velocity
    Particle(i).Velocity =zeros(VarSize);
    
    %checking cost function value
    Particle(i).Cost = CostFunction(Particle(i).Position,data_all);
    
    %update personal best
    Particle(i).Best.Position =Particle(i).Position;
    Particle(i).Best.Cost =Particle(i).Cost;
    
    %Update global best
    
    if Particle(i).Best.Cost < GlobalBest.Cost
        GlobalBest= Particle(i).Best;
    end
    
end

%Best cost value in each iterations
BestCost=zeros(MaxIt,1);

%%  Main loop of PSO

for it=1: MaxIt
    
    for i=1:nPoP
        
        %Update Velocity
        
        Particle(i).Velocity= w*Particle(i).Velocity+ ...
            c1*rand(VarSize).*(Particle(i).Best.Position-Particle(i).Position) ...
            + c2*rand(VarSize).*(GlobalBest.Position-Particle(i).Position);
        
        %Update Position
        
        Particle(i).Position=Particle(i).Position+Particle(i).Velocity;
        
        %cost for this iteration
        
        Particle(i).Cost = CostFunction(Particle(i).Position,data_all);
        
        %Update Personal Best
        
        if Particle(i).Cost < Particle(i).Best.Cost
            
            Particle(i).Best.Position =Particle(i).Position;
            Particle(i).Best.Cost =Particle(i).Cost;
            
            %Update Global Best
            
            if Particle(i).Best.Cost < GlobalBest.Cost
                GlobalBest= Particle(i).Best;
            end
            
        end
            
    end
    
    %fprintf('\t best density %f\n',GlobalBest.Position(end))
    
    %Display iteration information
    
    
    %w=w-((w-0.1)/100)*it;
    w=w*.99;
    %w=.9-(0.8*it)/MaxIt
    %Store Best Cost value
    BestCost(it)=myCostFunction(GlobalBest.Position,data_all);
    %pause(0.1)
    if BestCost(it)<0.005
        break
    end
                    
end
    data_value(cnt,:)= (GlobalBest.Position)';  

end
filename=sprintf('real_profile1_point%d_locations.dat',plg);
save(filename,'data_value','-ascii')
clear data_value
end
%%

function val=myCostFunction(x,data)
    
    sz=length(x);
    sz=sz-1;
    xx1=x(1,1:sz/2);
    yy1=x(1,(sz/2)+1:sz);
    
    %All observation points (here 0 to 2 km) of 40 data point
    x_obs=linspace(-6110.11236731062,6836.48208131634,82);
    x_obs=x_obs.*1;      %points in meter
    %z observation point at 0
    z_obs=0;
    
    [x1,y1]=poly_points(xx1,yy1);
       
    zz1=polygrav_arctan(x_obs,z_obs,x1,y1,x(sz+1));
    zz1=zz1*10^5;
    zz2=polygrad_zx(x_obs,z_obs,x1,y1,x(sz+1));
    zz2=zz2*10^8;
    d_len=length(data);
    %val1=norm(zz1-data(1,1:d_len/2));
    %val2=norm(zz2-data(1,(d_len/2)+1:d_len));
    zz=[zz1,zz2];
    val=norm(zz-data);
    %val=(100/(length(x_obs)))*sqrt(sum(((zz-data)./data).^2));
    %val=norm(zz-data)^2;

end

%function val=constrained1(x)

%    sz=length(x);
%    x1=x(1,1:sz/2);
%    y1=x(1,(sz/2)+1:sz);
    
%    perim=poly_perim(x1,y1);
%    gg=(perim-1500);
%    val=(max(0,gg))^2;
%end
function val=constrained1(x)

    sz=length(x);
    sz=sz-1;
    x1=x(1,1:sz/2);
    y1=x(1,(sz/2)+1:sz);
   
    gg=((poly_perim(x1,y1))^2/(4*pi*poly_area(x1,y1)))-1.4;
    val=(max(0,gg))^2;
end

function val=constrained2(x)

    sz=length(x);
    sz=sz-1;
    x1=x(1,1:sz/2);
    y1=x(1,(sz/2)+1:sz);
   
    gg=-(((poly_perim(x1,y1))^2/(4*pi*poly_area(x1,y1)))-1);
    val=(max(0,gg))^2;
end

function val=constrained3(x)

    sz=length(x);
    sz=sz-1;
    x1=x(1,1:sz/2);
    y1=x(1,(sz/2)+1:sz);
    
    area=poly_area(x1,y1);
    gg=(-area+350*350);
    val=(max(0,gg))^2;
end

function val=constrained4(x)

    sz=length(x);
    sz=sz-1;
    x1=x(1,1:sz/2);
    y1=x(1,(sz/2)+1:sz);
    
    area=poly_area(x1,y1);
    gg=(area-1000*1000);
    val=(max(0,gg))^2;
end

function val=constrained5(x)
    
    sz=length(x);
    sz=sz-1;
    xx1=x(1,1:sz/2);
    yy1=x(1,(sz/2)+1:sz);
    
    val1=0; val2=0;
    for i=1:sz/2
        gg1=-xx1(1,i)-2000;
        val1=val1+(max(0,gg1))^2;
        
        gg2=-yy1(1,i)-0;
        val2=val2+(max(0,gg2))^2;
    end
    val=val1+val2;
end

function val=constrained6(x)
    
    sz=length(x);
    sz=sz-1;
    xx1=x(1,1:sz/2);
    yy1=x(1,(sz/2)+1:sz);
    
    val1=0; val2=0;
    for i=1:sz/2
        %gg1=xx1(1,i)-2000;
        gg1=xx1(1,i)-2000;
        val1=val1+(max(0,gg1))^2;
        
        %gg2=yy1(1,i)-1000;
        gg2=yy1(1,i)-5000;
        val2=val2+(max(0,gg2))^2;
    end
    val=val1+val2;
end

function val=constrained7(x)
    rr=x(end);
    gg=(-rr-5000);
    val=(max(0,gg))^2;
end

function val=constrained8(x)
    rr=x(end);
    gg=(rr+2000);
    val=(max(0,gg))^2;
end