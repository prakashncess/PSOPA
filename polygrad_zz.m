function grav=polygrad_zz(x_obs,z_obs,x,z,roh)

%Polygrav function calculates z component of gravity field for any polygon 
%shape 2d body having finite density contrast. This program based on line
%integral in anticlockwise direction using Gauss Legendre quadrature
%integral formula. For more detail go through Zhou 2008. 

%x_obs is a vector containg observation points in x direction.
%z_obs level at which we are calculating gravity field, positive along
%downward direction
%x is the x coordinates of polygon body in counterclockwise direction. 
%z is the z coordinates of polygon body in counterclockwise direction. 
%roh is the density contrast.

% Always keep in mind x & z should always be taken in counter clockwise
% direction, otherwise sign convention will create problem while running. 

    n_poly=length(x); %length of the polygon
    x(length(x)+1)=x(1);% end point should be 1st point to close the integral
    z(length(z)+1)=z(1);% end point should be 1st point to close the integral
    G=6.67408*10^-11;% Gravitational constant in S.I.
    [t,c]=lgwt(10,0,1); %Legendre Gaussian quadrature integral points.
    for i=1:length(x_obs) % Loop for all observation points. 
        for j=1:n_poly % Loop for line integral over all sides of polygon 
            % Refer to Zhou 2008 paper for below steps, basically line
            % integral procedures. 
            
            ax1=(x(j).*(1-t)+x(j+1).*t-x_obs(i));
            ax2=(z(j).*(1-t)+z(j+1).*t-z_obs);
            ax=2.*roh.*G.*(ax1./(ax1.^2+ax2.^2)).*(z(j+1)-z(j));
            value(j) = sum(c.*ax);
            
        end
        grav(i)=sum(value(:)); % Combined gravity field for all points. 
    end
end

