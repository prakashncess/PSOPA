%Matlab code for plotting data 
clear all
close all
        dataa=importdata('model2_point6_without_noise_varying_rho_grd.dat');
        model=importdata('Model2.dat');
        noise=0; %0 for no noise 1 for noisy data
        sz_dt=size(dataa);
        tt=mod(sz_dt(2),2);
        nn=30; %how many best pso have to take  
        pt_sl=2;
                grd_wt=10;
                x_actual=model(1,:);
                y_actual=model(2,:);
                [x,y]=poly_points(x_actual,y_actual);
                %All observation points (here 0 to 2 km) of 40 data point
                x_obs=linspace(0,2,100);
                x_obs=x_obs.*1000;      %points in meter
                %z observation point at 0
                z_obs=0;
                rho=2700;
                zz1=polygrav_arctan(x_obs,z_obs,x,y,rho);
                data1=zz1*10^5;

                zz2=polygrad_zz(x_obs,z_obs,x,y,rho);
                data2=(zz2*10^9)/grd_wt;
                
                zz3=polygrad_zx(x_obs,z_obs,x,y,rho);
                data3=(zz3*10^9)/grd_wt;
                
                if noise==1
                    %adding of 10% Gaussian noise
                    %data1= data1+randn(size(data1)).*0.1.*data1;
                    %data2= data2+randn(size(data2)).*0.1.*data2;
                    %data3= data3+randn(size(data3)).*0.1.*data3;
                    noisedata1 = 0.1 * data1;
                    noise = noisedata1 .* randn(1, length(data3));
                    data1 = data1 + noise;
                    
                     noisedata2 = 0.1 * data2;
                    noise = noisedata2 .* randn(1, length(data2));
                    data2 = data2 + noise;
                    
                     noisedata3 = 0.1 * data3;
                    noise = noisedata1 .* randn(1, length(data3));
                    data3 = data3 + noise;
                end
                    data=[data1,data2,data3];
         if tt==0       
                for cnt=1:sz_dt(1)
                    xx_val=dataa(cnt,:);
                    sz=length(xx_val);
                    %sz=sz-1;
                    xx1=xx_val(1,1:sz/2);
                    yy1=xx_val(1,(sz/2)+1:sz);
                    [xt1,yt1]=poly_points(xx1,yy1);
                    [cmn_area,excd_area,xc_actual,yc_actual,xc_estmd,yc_estmd]=common_area(x_actual,y_actual,xx1,yy1);
                    commn_ara(cnt)=cmn_area; excdd_ara(cnt)=excd_area; xc_estm(cnt)=xc_estmd; yc_estm(cnt)=yc_estmd; 
                    gg_bst(cnt,:)=xx_val; %msft(cnt)=norm(data1-10^5*polygrav_arctan(x_obs,z_obs,xt1,yt1,rho));
                    msft(cnt)=100*sqrt((1/100).*sum((data1-10^5*polygrav_arctan(x_obs,z_obs,xt1,yt1,rho)).^2))./(max(data1(:))-min(data1(:)));
                end

               
            [ss,ii]=sort(commn_ara,'descend');

            xx_val=gg_bst(ii(pt_sl),:);
            xx1=xx_val(1,1:sz/2);
            yy1=xx_val(1,(sz/2)+1:sz);  
            misfit=msft(ii(pt_sl));
            [cmn_area,excd_area,xc_actual,yc_actual,xc_estmd,yc_estmd]=common_area(x_actual,y_actual,xx1,yy1);
            fprintf('For best estimation: Best Cost = %f\n',misfit)
            fprintf('\ncmn_area=%2.2f \texcc_area=%2.2f\n\t xc_actual=%2.2f yc_actual=%2.2f\n\t xc_estmd=%2.2f yc_estmd=%2.2f \n',cmn_area,excd_area,xc_actual,yc_actual,xc_estmd,yc_estmd);

            figure(1)
            clf
            subplot(4,1,4)
            hold on
            poly_plot(xx_val)
            fill(x,y,'b','facealpha',.25)
            xlim([0 2000])
            ylim([0 1000])
            ylabel('Depth (m)')
            xlabel('Distance (m)')
            title('Actual and inverted 2d Shape')
            %legend('Actual','Inverted','location','best')
            box on

            [x1,y1]=poly_points(xx1,yy1);
            zz1=polygrav_arctan(x_obs,z_obs,x1,y1,rho);
            zz1=zz1*10^5;

            zz2=polygrad_zz(x_obs,z_obs,x1,y1,rho);
            zz2=zz2*10^9;
            
            zz3=polygrad_zx(x_obs,z_obs,x1,y1,rho);
            zz3=zz3*10^9;

                    subplot(4,1,1)
                    hold on
                    plot(x_obs,data1,'b.')
                    plot(x_obs,zz1,'r')
                    xlabel('Distance (m)')
                    ylabel({'Gravity Anomaly'; '(mGal)'})
                    title('Observed and inverted field')
                    legend('Obsrved','Inverted','location','best')
                    box on

                    subplot(4,1,2)
                    hold on
                    plot(x_obs,data2*grd_wt,'b.')
                    plot(x_obs,zz2,'r')
                    xlabel('Distance (m)')
                    ylabel({'Gradient G_{zz}'; '(Eotvos)'})
                    title('Observed and inverted field')
                    legend('Obsrved','Inverted','location','best')
                    box on
                    
                    subplot(4,1,3)
                    hold on
                    plot(x_obs,data3*grd_wt,'b.')
                    plot(x_obs,zz3,'r')
                    xlabel('Distance (m)')
                    ylabel({'Gradient G_{zx}'; '(Eotvos)'})
                    title('Observed and inverted field')
                    legend('Obsrved','Inverted','location','best')
                    box on
        %misfit=norm(data1-zz1)+norm(data2*grd_wt-zz2)+norm(data3*grd_wt-zz3);
        commn_ara=commn_ara(ii(1:nn));
        excdd_ara=excdd_ara(ii(1:nn));
        xc_estm=xc_estm(ii(1:nn));
        yc_estm=yc_estm(ii(1:nn));
        msft=msft(ii(1:nn));
        
        mn_cmm_ara= mean(commn_ara); std_cmm_ara=std(commn_ara);
        fprintf('\tMean common area=%f \t Standard Deviation=%2.2f\n',mn_cmm_ara,std_cmm_ara)
        mn_exs_ara= mean(excdd_ara); std_exs_ara=std(excdd_ara);
        fprintf('\tMean excess area=%f \t Standard Deviation=%2.2f\n',mn_exs_ara,std_exs_ara)
        mn_xc=mean(xc_estm);  std_xc=std(xc_estm);
        fprintf('\tMean x position=%f \t Standard Deviation=%2.2f\n',mn_xc,std_xc)
        mn_yc=mean(yc_estm);  std_yc=std(yc_estm);
        fprintf('\tMean z position=%f \t Standard Deviation=%2.2f\n',mn_yc,std_yc)
        mn_msft=mean(msft);  std_msft=std(msft);
        fprintf('\tMean misfit=%f \t Standard Deviation=%2.2f\n',mn_msft,std_msft)
        
        edge_area=linspace(min(commn_ara),max(commn_ara),10);
        edge_xarea=linspace(min(excdd_ara),max(excdd_ara),10);
        edge_xc=linspace(min(xc_estm),max(xc_estm),10);
        edge_yc=linspace(min(yc_estm),max(yc_estm),10);
        N_area = histcounts(commn_ara',edge_area);
        N_xarea = histcounts(excdd_ara',edge_xarea);
        N_xc = histcounts(xc_estm',edge_xc);
        N_yc = histcounts(yc_estm',edge_yc);



        figure(2)
        subplot(1,4,1)
            binWidth = edge_area(2)-edge_area(1);
            binCenters = edge_area(1:end-1) + binWidth/2;
            bar(binCenters, (N_area/sum(N_area))*100, 'BarWidth', 0.5, 'FaceColor', [0 .5 .5],'EdgeColor','k','LineWidth',1.3);
            grid on;
            xlabel('Overlaping area (m^2)')
            ylabel('Frequencies')
        subplot(1,4,2)
            binWidth = edge_xarea(2)-edge_xarea(1);
            binCenters = edge_xarea(1:end-1) + binWidth/2;
            bar(binCenters, (N_xarea/sum(N_xarea))*100, 'BarWidth', 0.5, 'FaceColor', [0 .5 .5],'EdgeColor','k','LineWidth',1.3);
            grid on;
            xlabel('Excess area (m^2)')
            ylabel('Frequencies') 
            %title(sprintf('Frequency of parameters value for %d point polygon',plg))
        subplot(1,4,3)    
            binWidth2 = edge_xc(2)-edge_xc(1);
            binCenters2 = edge_xc(1:end-1) + binWidth2/2;
            bar(binCenters2, (N_xc/sum(N_xc))*100, 'BarWidth', 0.5, 'FaceColor', [0 .5 .5],'EdgeColor','k','LineWidth',1.3);
            grid on;
            xlabel('Location (m)')
            ylabel('Frequencies')

        subplot(1,4,4)  
            binWidth = edge_yc(2)-edge_yc(1);
            binCenters = edge_yc(1:end-1) + binWidth/2;
            bar(binCenters, (N_yc/sum(N_yc))*100, 'BarWidth', 0.5, 'FaceColor', [0 .5 .5],'EdgeColor','k','LineWidth',1.3);
            grid on;
            xlabel('Depth (m)')
            ylabel('Frequencies')
            
        else
                
                for cnt=1:sz_dt(1)
                    xx_val=dataa(cnt,:);
                    sz=length(xx_val);
                    sz=sz-1;
                    xx1=xx_val(1,1:sz/2);
                    yy1=xx_val(1,(sz/2)+1:sz);
                    [xt1,yt1]=poly_points(xx1,yy1);
                    [cmn_area,excd_area,xc_actual,yc_actual,xc_estmd,yc_estmd]=common_area(x_actual,y_actual,xx1,yy1);
                    commn_ara(cnt)=cmn_area; excdd_ara(cnt)=excd_area; xc_estm(cnt)=xc_estmd; yc_estm(cnt)=yc_estmd; 
                    rho_estd(cnt)=xx_val(end);gg_bst(cnt,:)=xx_val; 
                    %msft(cnt)=norm(10^5*polygrav_arctan(x_obs,z_obs,x,y,rho)-10^5*polygrav_arctan(x_obs,z_obs,xt1,yt1,rho_estd(cnt)));
                    msft(cnt)=100*(sqrt((1/100).*sum((data1-10^5*polygrav_arctan(x_obs,z_obs,xt1,yt1,rho_estd(cnt))).^2)))./(max(data1(:))-min(data1(:)));
                end

            [ss,ii]=sort(commn_ara,'descend');
            xx_val=gg_bst(ii(pt_sl),:);
            xy_val=gg_bst(ii(pt_sl),1:end-1);
            xx1=xy_val(1,1:sz/2);
            yy1=xy_val(1,(sz/2)+1:sz);  
            rho=rho_estd(ii(pt_sl));
            misfit=msft(ii(pt_sl));
            [cmn_area,excd_area,xc_actual,yc_actual,xc_estmd,yc_estmd]=common_area(x_actual,y_actual,xx1,yy1);
            fprintf('For best estimation: Best Cost = %f\n',misfit)
            fprintf('\ncmn_area=%2.2f \texcc_area=%2.2f\n\t xc_actual=%2.2f yc_actual=%2.2f\n\t xc_estmd=%2.2f yc_estmd=%2.2f\n\t rho_actual=2700 rho_estmtd=%f \n\n',cmn_area,excd_area,xc_actual,yc_actual,xc_estmd,yc_estmd,rho);

            figure(1)
            clf
            subplot(4,1,4)
            hold on
            poly_plot(xy_val)
            fill(x,y,'b','facealpha',.25)
            xlim([0 2000])
            ylim([0 1000])
            ylabel('Depth (m)')
            xlabel('Distance (m)')
            title('Actual and inverted 2d Shape')
            %legend('Actual','Inverted','location','best')
            box on

            [x1,y1]=poly_points(xx1,yy1);
            zz1=polygrav_arctan(x_obs,z_obs,x1,y1,rho);
            zz1=zz1*10^5;

            zz2=polygrad_zz(x_obs,z_obs,x1,y1,rho);
            zz2=zz2*10^9;
            
            zz3=polygrad_zx(x_obs,z_obs,x1,y1,rho);
            zz3=zz3*10^9;

                    subplot(4,1,1)
                    hold on
                    plot(x_obs,data1,'b.')
                    plot(x_obs,zz1,'r')
                    xlabel('Distance (m)')
                    ylabel({'Gravity Anomaly'; '(mGal)'})
                    title('Observed and inverted field')
                    legend('Obsrved','Inverted','location','best')
                    box on

                    subplot(4,1,2)
                    hold on
                    plot(x_obs,data2*grd_wt,'b.')
                    plot(x_obs,zz2,'r')
                    xlabel('Distance (m)')
                    ylabel({'Gradient G_{zz}'; '(Eotvos)'})
                    title('Observed and inverted field')
                    legend('Obsrved','Inverted','location','best')
                    box on
                    
                    subplot(4,1,3)
                    hold on
                    plot(x_obs,data3*grd_wt,'b.')
                    plot(x_obs,zz3,'r')
                    xlabel('Distance (m)')
                    ylabel({'Gradient G_{zx}'; '(Eotvos)'})
                    title('Observed and inverted field')
                    legend('Obsrved','Inverted','location','best')
                    box on
                    
                      
            commn_ara=commn_ara(ii(1:nn));
            excdd_ara=excdd_ara(ii(1:nn));
            rho_estd=rho_estd(ii(1:nn));
            xc_estm=xc_estm(ii(1:nn));
            yc_estm=yc_estm(ii(1:nn));
            msft=msft(ii(1:nn));
            mn_cmm_ara= mean(commn_ara); std_cmm_ara=std(commn_ara);
            fprintf('\tMean common area=%f \t Standard Deviation=%2.2f\n',mn_cmm_ara,std_cmm_ara)
            mn_exs_ara= mean(excdd_ara); std_exs_ara=std(excdd_ara);
            fprintf('\tMean excess area=%f \t Standard Deviation=%2.2f\n',mn_exs_ara,std_exs_ara)
            mn_rho_estd= mean(rho_estd); std_rho_estd=std(rho_estd);
            fprintf('\tMean density=%f \t Standard Deviation=%2.2f\n',mn_rho_estd,std_rho_estd)
            mn_xc=mean(xc_estm);  std_xc=std(xc_estm);
            fprintf('\tMean x position=%f \t Standard Deviation=%2.2f\n',mn_xc,std_xc)
            mn_yc=mean(yc_estm);  std_yc=std(yc_estm);
            fprintf('\tMean z position=%f \t Standard Deviation=%2.2f\n',mn_yc,std_yc)
            mn_msft=mean(msft);  std_msft=std(msft);
            fprintf('\tMean misfit=%f \t Standard Deviation=%2.2f\n',mn_msft,std_msft)
            edge_area=linspace(min(commn_ara),max(commn_ara),10);
            edge_xarea=linspace(min(excdd_ara),max(excdd_ara),10);
            edge_rho=linspace(min(rho_estd),max(rho_estd),10);
            edge_xc=linspace(min(xc_estm),max(xc_estm),10);
            edge_yc=linspace(min(yc_estm),max(yc_estm),10);
            N_area = histcounts(commn_ara',edge_area);
            N_xarea = histcounts(excdd_ara',edge_xarea);
            N_rho= histcounts(rho_estd',edge_rho);
            N_xc = histcounts(xc_estm',edge_xc);
            N_yc = histcounts(yc_estm',edge_yc);
            figure(2)
            %title(sprintf('Frequency of parameters value for %d point polygon',plg))        
            subplot(1,5,1)
                binWidth = edge_area(2)-edge_area(1);
                binCenters = edge_area(1:end-1) + binWidth/2;
                bar(binCenters, (N_area/sum(N_area))*100, 'BarWidth', 0.5, 'FaceColor', [0 .5 .5],'EdgeColor','k','LineWidth',1.3);
                grid on;
                xlabel({'Overlaping area'; '(m^2)'})
                ylabel('Frequencies')

            subplot(1,5,2)
                binWidth = edge_xarea(2)-edge_xarea(1);
                binCenters = edge_xarea(1:end-1) + binWidth/2;
                bar(binCenters, (N_xarea/sum(N_xarea))*100, 'BarWidth', 0.5, 'FaceColor', [0 .5 .5],'EdgeColor','k','LineWidth',1.3);
                grid on;
                xlabel({'Excess area';'(m^2)'})
                ylabel('Frequencies')    
            subplot(1,5,3)    
                binWidth2 = edge_xc(2)-edge_xc(1);
                binCenters2 = edge_xc(1:end-1) + binWidth2/2;
                bar(binCenters2, (N_xc/sum(N_xc))*100, 'BarWidth', 0.5, 'FaceColor', [0 .5 .5],'EdgeColor','k','LineWidth',1.3);
                grid on;
                xlabel('Location (m)')
                ylabel('Frequencies')
                %title(sprintf('Frequency of parameters value for %d point polygon',plg))  
            subplot(1,5,4)  
                binWidth = edge_yc(2)-edge_yc(1);
                binCenters = edge_yc(1:end-1) + binWidth/2;
                bar(binCenters, (N_yc/sum(N_yc))*100, 'BarWidth', 0.5, 'FaceColor', [0 .5 .5],'EdgeColor','k','LineWidth',1.3);
                grid on;
                xlabel('Depth (m)')
                ylabel('Frequencies')
             subplot(1,5,5)  
                binWidth = edge_rho(2)-edge_rho(1);
                binCenters = edge_rho(1:end-1) + binWidth/2;
                bar(binCenters, (N_rho/sum(N_rho))*100, 'BarWidth', 0.5, 'FaceColor', [0 .5 .5],'EdgeColor','k','LineWidth',1.3);
                grid on;
                xlabel('Density (kg/m^3)')
                ylabel('Frequencies')
        end