%Matlab function for finding common area percentage and center positions of
%two polygons

function [cmn_area,excd_area,xc_actual,yc_actual,xc_estmd,yc_estmd]=common_area(x_actual,y_actual,x_estmd,y_estmd)

    %Outputs     cmn_area= percentage of common intesected area
    %            excd_area=percentage of excess area outside the actual
    %            polygon
    %           (xc_actual,yc_actual)= Center position of actual polygon
    %           (xc_estmd,yc_estmd)  = Center positon of inverted polygon
    
    %Inputs     x_actual = x positions of the actual polygon
    %           y_actual = y positions of the actual polygon
    %           x_estmd = x positions of the estimated polygon
    %           y_estmd = y positions of the estimated polygon
    
                % location of inverted polygon
                xx1=x_estmd;
                yy1=y_estmd;
                
                %location of actual polygon
                xx2=x_actual;
                yy2=y_actual;

                [x1,y1]=poly_points(xx1,yy1);
                [x2,y2]=poly_points(xx2,yy2);
                
                
                poly1= polyshape(x1,y1);
                poly2= polyshape(x2,y2);
                
                %location of intersection points of both polygons
                polyout = intersect(poly1,poly2);
                
                %Vertices of intersecting points 
                xy_intr=polyout.Vertices;

                %center position of actual and estimated polygon                
                [xc_estmd,yc_estmd] = centroid(poly1);
                [xc_actual,yc_actual] = centroid(poly2);
                
                xcc=squeeze(xy_intr(:,1)); ycc=squeeze(xy_intr(:,2));
                tf=isempty(xcc);

                if tf==1
                    comon_area=0;
                else
                    comon_area=poly_area(xcc',ycc');
                end

                
                actual_area=poly_area(xx2,yy2);
                %estimated_area=poly_area(xx1,yy1);
                %outside_area=abs(common_area-estimated_area);

                cmn_area=(comon_area/actual_area)*100;
                excd_area=((poly_area(x_estmd,y_estmd)-comon_area)/actual_area)*100;
                
end

    