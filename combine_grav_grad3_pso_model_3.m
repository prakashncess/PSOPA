%%Matlab code for Particle Swarm Algorithm 
clear all
close all
%Three different types of Model
%model for creating sphere
 rad=250; y_cen=500; x_cen=1000;
 tht=linspace(0,2*pi,20);
 plg=4; %total number of vertex of the polygon    
 %creting bigger rectangle having same dim
 x_actual=rad*cos(tht(1:end-1))+x_cen;
 y_actual=rad*sin(tht(1:end-1))+y_cen;
model2=[x_actual;y_actual];
save('Model2.dat','model2','-ascii')

for ww=1:2          %ww=1 for without noise or ww=2 for with noise
        for plg=4:6 %total number of vertex of the polygon
            
            for rhh=1:2% rhh=1 for fixed rho rhh=2 for varying rho
                %total 65 generations have been considered
                for cnt=1:65
                    fprintf('ww=%d , rhh= %d, numberof polygon=%d pso_count=%d\n',ww,rhh,plg,cnt)
                    [x,y]=poly_points(x_actual,y_actual);
                    %gravity field for this bigger domain
                    %All observation points (here 0 to 2 km) of 40 data point
                    x_obs=linspace(0,2,100);
                    x_obs=x_obs.*1000;      %points in meter
                    %z observation point at 0
                    z_obs=0;
                    %density in kg/m^3
                    rho=2700;
                    %gravity field 
                    zz1=polygrav_arctan(x_obs,z_obs,x,y,rho);
                    data1=zz1*10^5;
                    %horizontal gradient
                    zz2=polygrad_zz(x_obs,z_obs,x,y,rho);
                    data2=zz2*10^8;
                    zz3=polygrad_zx(x_obs,z_obs,x,y,rho);
                    data3=zz3*10^8;
                    %INCORPORATING NOISE
                    if ww==2
                        %adding of 10% Gaussian noise
                        data1= data1+randn(size(data1)).*0.1.*data1;
                        data2= data2+randn(size(data2)).*0.1.*data2;
                        data3= data3+randn(size(data3)).*0.1.*data3;
                    end
                    data=[data1,data2,data3];
                    if rhh==1
       %Starting of PSO algorithm 
                        %%Problem Definition
                        CostFunction1 =@(x,data,rho) myCostFunction1(x,data,rho)+...
                            1000*(constrained11(x)+constrained12(x)+constrained13(x)+constrained14(x)+constrained15(x)+constrained16(x));

                        %nVar=8;                    %Number of Unknown Variable 
                        nVar=plg*2;                 %Number of Unknown Variable 

                        VarSize = [1 nVar];     %Matrix size of Decision variables

                        VarMin= -ones(1,nVar);          %Lower Bound of Unknown variable
                        VarMax= ones(1,nVar);           %Upper Bound of Unknown variable

                        %% Parameters of PSO

                        MaxIt = 300;            %Maximum number of iterations

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
                            Particle(i).Cost = CostFunction1(Particle(i).Position,data,rho);

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

                                Particle(i).Cost = CostFunction1(Particle(i).Position,data,rho);

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

                            %Store Best Cost value
                            BestCost(it)=GlobalBest.Cost;

                            %Display iteration information
                            %disp (['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))])


                            %w=w-((w-0.1)/100)*it;
                            w=w*.99;
                            %w=.9-(0.8*it)/MaxIt

                            BestCost(it)=myCostFunction1(GlobalBest.Position,data,rho);
                            %pause(0.1)
                            if BestCost(it)<0.005
                                break
                            end

                        end
                    else
                        CostFunction2 =@(x,data) myCostFunction2(x,data)+...
                        1000*(constrained21(x)+constrained22(x)+constrained23(x)+constrained24(x)+constrained25(x)+constrained26(x)+...
                              constrained27(x)+constrained28(x));

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
                            Particle(i).Cost = CostFunction2(Particle(i).Position,data);

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

                                Particle(i).Cost = CostFunction2(Particle(i).Position,data);

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

                            %Store Best Cost value
                            BestCost(it)=myCostFunction2(GlobalBest.Position,data);
                            %fprintf('rho=%f\n',GlobalBest.Position(end))
                            %Display iteration information


                            %w=w-((w-0.1)/100)*it;
                            w=w*.99;
                            %w=.9-(0.8*it)/MaxIt

                            %pause(0.1)
                            if GlobalBest.Cost<0.005
                                break
                            end
                        end
                    end

                    data_value(cnt,:)= (GlobalBest.Position)';  

                end
                if ww==1
                   if rhh==1
                        filename=sprintf('model2_point%d_without_noise_fixed_rho_grd.dat',plg);
                    else
                        filename=sprintf('model2_point%d_without_noise_varying_rho_grd.dat',plg);
                    end
                else
                    if rhh==1
                        filename=sprintf('model2_point%d_with_noise_fixed_rho_grd.dat',plg);
                    else
                        filename=sprintf('model2_point%d_with_noise_varying_rho_grd.dat',plg);
                    end
                end
                    save(filename,'data_value','-ascii')
                    clear data_value
            end
        end

end
%%    
function val=myCostFunction1(x,data,rho)
    
    sz=length(x);
    %sz=sz-1;
    xx1=x(1,1:sz/2);
    yy1=x(1,(sz/2)+1:sz);
    
    %All observation points (here 0 to 2 km) of 40 data point
    x_obs=linspace(0,2,100);
    x_obs=x_obs.*1000;      %points in meter
    %z observation point at 0
    z_obs=0;
    
    [x1,y1]=poly_points(xx1,yy1);
      
    zz1=polygrav_arctan(x_obs,z_obs,x1,y1,rho);
    zz1=zz1*10^5;
    zz2=polygrad_zz(x_obs,z_obs,x1,y1,rho);
    zz2=zz2*10^8;
    zz3=polygrad_zx(x_obs,z_obs,x1,y1,rho);
    zz3=zz3*10^8;
    zz=[zz1,zz2,zz3];
    val=norm(zz-data);
end

function val=myCostFunction2(x,data)
    
    sz=length(x);
    sz=sz-1;
    xx1=x(1,1:sz/2);
    yy1=x(1,(sz/2)+1:sz);
    
    %All observation points (here 0 to 2 km) of 40 data point
    x_obs=linspace(0,2,100);
    x_obs=x_obs.*1000;      %points in meter
    %z observation point at 0
    z_obs=0;
    
    [x1,y1]=poly_points(xx1,yy1);
       
    zz1=polygrav_arctan(x_obs,z_obs,x1,y1,x(sz+1));
    zz1=zz1*10^5;
    zz2=polygrad_zz(x_obs,z_obs,x1,y1,x(sz+1));
    zz2=zz2*10^8;
    zz3=polygrad_zx(x_obs,z_obs,x1,y1,x(sz+1));
    zz3=zz3*10^8;
    zz=[zz1,zz2,zz3];
    val=norm(zz-data);

end

%function val=constrained1(x)

%    sz=length(x);
%    x1=x(1,1:sz/2);
%    y1=x(1,(sz/2)+1:sz);
    
%    perim=poly_perim(x1,y1);
%    gg=(perim-1500);
%    val=(max(0,gg))^2;
%end
function val=constrained11(x)

    sz=length(x);
    %sz=sz-1;
    x1=x(1,1:sz/2);
    y1=x(1,(sz/2)+1:sz);
   
    gg=((poly_perim(x1,y1))^2/(4*pi*poly_area(x1,y1)))-1.3;
    val=(max(0,gg))^2;
end

function val=constrained12(x)

    sz=length(x);
    %sz=sz-1;
    x1=x(1,1:sz/2);
    y1=x(1,(sz/2)+1:sz);
   
    gg=-(((poly_perim(x1,y1))^2/(4*pi*poly_area(x1,y1)))-1);
    val=(max(0,gg))^2;
end

function val=constrained13(x)

    sz=length(x);
    %sz=sz-1;
    x1=x(1,1:sz/2);
    y1=x(1,(sz/2)+1:sz);
    
    area=poly_area(x1,y1);
    gg=(-area+2*2);
    val=(max(0,gg))^2;
end

function val=constrained14(x)

    sz=length(x);
    %sz=sz-1;
    x1=x(1,1:sz/2);
    y1=x(1,(sz/2)+1:sz);
    
    area=poly_area(x1,y1);
    gg=(area-500*500);
    val=(max(0,gg))^2;
end

function val=constrained15(x)
    
    sz=length(x);
    %sz=sz-1;
    xx1=x(1,1:sz/2);
    yy1=x(1,(sz/2)+1:sz);
    
    val1=0; val2=0;
    for i=1:sz/2
        gg1=-xx1(1,i)+100;
        val1=val1+(max(0,gg1))^2;
        
        gg2=-yy1(1,i)+100;
        val2=val2+(max(0,gg2))^2;
    end
    val=val1+val2;
end

function val=constrained16(x)
    
    sz=length(x);
    %sz=sz-1;
    xx1=x(1,1:sz/2);
    yy1=x(1,(sz/2)+1:sz);
    
    val1=0; val2=0;
    for i=1:sz/2
        %gg1=xx1(1,i)-2000;
        gg1=xx1(1,i)-2000;
        val1=val1+(max(0,gg1))^2;
        
        %gg2=yy1(1,i)-1000;
        gg2=yy1(1,i)-750;
        val2=val2+(max(0,gg2))^2;
    end
    val=val1+val2;
end

function val=constrained21(x)

    sz=length(x);
    sz=sz-1;
    x1=x(1,1:sz/2);
    y1=x(1,(sz/2)+1:sz);
   
    gg=((poly_perim(x1,y1))^2/(4*pi*poly_area(x1,y1)))-1.3;
    val=(max(0,gg))^2;
end

function val=constrained22(x)

    sz=length(x);
    sz=sz-1;
    x1=x(1,1:sz/2);
    y1=x(1,(sz/2)+1:sz);
   
    gg=-(((poly_perim(x1,y1))^2/(4*pi*poly_area(x1,y1)))-1);
    val=(max(0,gg))^2;
end

function val=constrained23(x)

    sz=length(x);
    sz=sz-1;
    x1=x(1,1:sz/2);
    y1=x(1,(sz/2)+1:sz);
    
    area=poly_area(x1,y1);
    gg=(-area+2*2);
    val=(max(0,gg))^2;
end

function val=constrained24(x)

    sz=length(x);
    sz=sz-1;
    x1=x(1,1:sz/2);
    y1=x(1,(sz/2)+1:sz);
    
    area=poly_area(x1,y1);
    gg=(area-500*500);
    val=(max(0,gg))^2;
end

function val=constrained25(x)
    
    sz=length(x);
    sz=sz-1;
    xx1=x(1,1:sz/2);
    yy1=x(1,(sz/2)+1:sz);
    
    val1=0; val2=0;
    for i=1:sz/2
        gg1=-xx1(1,i)+100;
        val1=val1+(max(0,gg1))^2;
        
        gg2=-yy1(1,i)+100;
        val2=val2+(max(0,gg2))^2;
    end
    val=val1+val2;
end

function val=constrained26(x)
    
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
        gg2=yy1(1,i)-750;
        val2=val2+(max(0,gg2))^2;
    end
    val=val1+val2;
end

function val=constrained27(x)
    rr=x(end);
    gg=(rr-3000);
    val=(max(0,gg))^2;
end

function val=constrained28(x)
    rr=x(end);
    gg=-(rr-2400);
    val=(max(0,gg))^2;
end
          