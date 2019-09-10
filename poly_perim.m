%Matlab function for finding perimeter from randomly distributed points

function perim=poly_perim(x,y)

    %x is any points in x coordinates; 
    %y is any points in y coordinates;

    P=[x;y];
    c=mean(P,2);
    d=P-c;
    th=atan2(d(2,:),d(1,:));
    [th, idx]=sort(th);
    P=P(:,idx);
    P=[P P(:,1)];
    %fill(P(1,:),P(2,:),'.-r')
    %hold on
    %plot(x,y,'b*')

    perim=0;
    %finding area of the polygon from arranged random points
    for i=1:length(P)-1

        X1=P(1,i); X2=P(1,i+1);
        Y1=P(2,i); Y2=P(2,i+1);

        perim=perim+sqrt((X1-X2)^2+(Y1-Y2)^2);

    end
    
end