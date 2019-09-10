%Matlab code for plotting data 
clear all
close all
pt_sl=5; nn=40;
yy=importdata('C:\Users\Arka\Desktop\pso_ulti_final\real_data\Karrbo_Residual_Anomaly.bln');
x_obs1=yy(:,1);
data_grv=yy(:,2);
x_obs=linspace(min(x_obs1),max(x_obs1),length(yy));
data1=spline(x_obs1',data_grv',x_obs);
z_obs=0;
data1=data1*10^-5;
data2=diff(data1(:))./diff(x_obs(:));
data2=data2*10^9;
data3=spline(x_obs(1:end-1),data2,x_obs);
data_gv=data1*10^5;
data_gr=data3;
dataa=importdata('real_profile6_point7_locations.dat');
sz_dt=size(dataa);
for cnt=1:sz_dt(1)
    
    xx_val=dataa(cnt,:);
    sz=length(xx_val);
    sz=sz-1;
    xx1=xx_val(1,1:sz/2);
    yy1=xx_val(1,(sz/2)+1:sz);
    [xv1,yv1]=poly_points(xx1,yy1);
    poly1= polyshape(xv1,yv1);
    [xc_estmdd(cnt),yc_estmdd(cnt)] = centroid(poly1);
    rho_estd(cnt)=xx_val(end);gg_bst(cnt,:)=xx_val; 
    area(cnt)=poly_area(xv1,yv1);
    ms1=100*sqrt((1/length(yy)).*sum((data_grv'-10^5*polygrav_arctan(x_obs,z_obs,xv1,yv1,rho_estd(cnt))).^2))./(max(data_grv(:))-min(data_grv(:)));
    ms2=100*sqrt((1/length(yy)).*sum((data3-10^9*polygrad_zx(x_obs,z_obs,xv1,yv1,rho_estd(cnt))).^2))./(max(data3(:))-min(data3(:)));
    msft(cnt)=(ms1+ms2)/2;
    %msft(cnt)=norm(10^5*polygrav_arctan(x_obs,z_obs,x,y,rho)-10^5*polygrav_arctan(x_obs,z_obs,xt1,yt1,rho_estd(cnt)));
    %msft(cnt)=100*sqrt((1/length(yy)).*sum((data_grv'-10^5*polygrav_arctan(x_obs,z_obs,xv1,yv1,rho_estd(cnt))).^2))./(max(abs(data_grv(:)))-min(abs(data_grv(:))));
    %msft(cnt)=100*(1/length(yy)).*sqrt(sum((abs(data_grv')-abs(10^5*polygrav_arctan(x_obs,z_obs,xv1,yv1,rho_estd(cnt)))./abs(data_grv')).^2));
end
            [ss,ii]=sort(msft);
            xx_val=gg_bst(ii(pt_sl),:);
            xy_val=gg_bst(ii(pt_sl),1:end-1);
            xx1=xy_val(1,1:sz/2);
            yy1=xy_val(1,(sz/2)+1:sz);  
            rho=rho_estd(ii(pt_sl));
            misfit=msft(ii(pt_sl));
            arrea=area(ii(pt_sl));
            xc_estmd=xc_estmdd(ii(pt_sl));
            yc_estmd=yc_estmdd(ii(pt_sl));
            fprintf('For best estimation: Best Cost = %f\n',misfit)
            fprintf('\n\tarea=%2.2f \n\txc_estmd=%2.2f yc_estmd=%2.2f\n\t rho_estmtd=%f \n\n',arrea,xc_estmd,yc_estmd,rho);
            
            figure(1)
            clf
            subplot(3,1,3)
            hold on
            poly_plot(xy_val)
            xlim([-13 13])
            ylim([0 20])
            ylabel('Depth (m)')
            xlabel('Distance (m)')
            title('Inverted 2d Shape')
            %legend('Actual','Inverted','location','best')
            box on

            [x1,y1]=poly_points(xx1,yy1);
            zz1=polygrav_arctan(x_obs,z_obs,x1,y1,rho);
            zz1=zz1*10^5;

                      
            zz2=polygrad_zx(x_obs,z_obs,x1,y1,rho);
            zz2=zz2*10^9;

                    subplot(3,1,1)
                    hold on
                    plot(x_obs,data_gv,'b.')
                    plot(x_obs,zz1,'r')
                    xlabel('Distance (m)')
                    ylabel({'Gravity Anomaly'; '(mGal)'})
                    title('Observed and Inverted field')
                    legend('Obsrved','Inverted','location','best')
                    box on

                    subplot(3,1,2)
                    hold on
                    plot(x_obs,data3,'b.')
                    plot(x_obs,zz2,'r')
                    xlabel('Distance (m)')
                    ylabel({'Gradient G_{zx}'; '(Eotvos)'})
                    title('Observed and Inverted field')
                    legend('Obsrved','Inverted','location','best')
                    box on
                    
                    
                      
            area=area(ii(1:nn));
            rho_estd=rho_estd(ii(1:nn));
            xc_estmdd=xc_estmdd(ii(1:nn));
            yc_estmdd=yc_estmdd(ii(1:nn));
            msft=msft(ii(1:nn));
            mn_cmm_ara= mean(area); std_cmm_ara=std(area);
            fprintf('\tMean area=%f \t Standard Deviation=%2.2f\n',mn_cmm_ara,std_cmm_ara)
            mn_rho_estd= mean(rho_estd); std_rho_estd=std(rho_estd);
            fprintf('\tMean density=%f \t Standard Deviation=%2.2f\n',mn_rho_estd,std_rho_estd)
            mn_xc=mean(xc_estmdd);  std_xc=std(xc_estmdd);
            fprintf('\tMean x position=%f \t Standard Deviation=%2.2f\n',mn_xc,std_xc)
            mn_yc=mean(yc_estmdd);  std_yc=std(yc_estmdd);
            fprintf('\tMean z position=%f \t Standard Deviation=%2.2f\n',mn_yc,std_yc)
            mn_msft=mean(msft);  std_msft=std(msft);
            fprintf('\tMean misfit=%f \t Standard Deviation=%2.2f\n',mn_msft,std_msft)
            
            