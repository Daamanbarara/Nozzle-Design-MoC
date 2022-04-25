%% Method of Characteristics Code for 2D-Minimum Length Nozzle
% Author- D. Barara

clc;
clear all;
%% Inputs

M_e=2.2;        % Exit Mach Number
gamma=1.4;      % Ratio of specific heats
n=50;           % Number of C_ characteristics (Right running characteristics)
x0=0;
y0=1;

%% Initializing arrays

nodes=(2+(n+1))*n/2;
node_no=zeros(n,1); 
K_n=zeros(n,1);        % K+ constant for C+ characteristic line (left running) 
K_p=zeros(n,1);        % K- constant for C- characteristic line (right running)
theta_n=zeros(n,1);    % Flow dfelection angle
nu_n=zeros(n,1);       % Prandtl Meyer function
mu_n=zeros(n,1);       % Mach angle
M_n=zeros(n,1);        % Mach number
x_n=zeros(n,1);        % X-coordinate position
y_n=zeros(n,1);        % Y-coordinate position

%% Computing Maximum wall angle and initial angle

theta_max=Prandtl_Meyer(M_e,gamma)/2;
% theta_initial=1-abs(((theta_max)-round(theta_max)));
if abs(((theta_max)-round(theta_max)))>=0.5
    theta_initial=abs(((theta_max)-round(theta_max)));
elseif abs(((theta_max)-round(theta_max)))<0.5
    theta_initial=1-abs(((theta_max)-round(theta_max)));
end
theta_delta=(theta_max-theta_initial)/(n-1);

%% Computation of variables at different nodes on the Characteristic lines

for j=1:n
    for i=1:n-j+1
        if j==1   % First Characteristic line
            theta_n(i,j)=theta_initial + (i-1)*theta_delta;
            nu_n(i,j)= theta_n(i,j);
            M_n(i,j)=Inverse_Prandtl_Meyer(nu_n(i,j),gamma);
            mu_n(i,j)=asind(1/M_n(i,j));
            K_n(i,j)= theta_n(i,j)+nu_n(i,j);
            K_p(i,j)= theta_n(i,j)-nu_n(i,j);
        elseif j>1  % Rest of the Characteristic data
            K_n(i,j)= K_n(i+1,j-1);
            if i==1 % Centerline
                K_p(i,j)= -1*K_n(i,j);
                theta_n(i,j)= 0;
                nu_n(i,j)= K_n(i,j);
                M_n(i,j)=Inverse_Prandtl_Meyer(nu_n(i,j),gamma);
                mu_n(i,j)=asind(1/M_n(i,j));
            else    % Internal Points (i>1)
                K_p(i,j)= K_p(i-1,j);
                theta_n(i,j)=(K_n(i+1,j-1)+K_p(i-1,j))/2;
                nu_n(i,j)=(K_n(i+1,j-1)-K_p(i-1,j))/2;
                M_n(i,j)=Inverse_Prandtl_Meyer(nu_n(i,j),gamma);
                mu_n(i,j)=asind(1/M_n(i,j));
            end
        end
    end
end

%% Coordinates Computation
% Coordinates of the charateristic line
for j=1:n
    for i=1:n-j+1
        if j==1
            if i==1
                A=((theta_n(i,j)))-((mu_n(i,j)));
                B=0;
                x_n(i,j)=(y0-tand(A)*(x0))/(tand(B)-tand(A));
                y_n(i,j)=0;
            else % (i>1) Rest of the points of First Characteristic line
                A=theta_n(i,j)-mu_n(i,j);
                B=((theta_n(i,j)+theta_n(i-1,j))/2)+((mu_n(i,j)+mu_n(i-1,j))/2);
                x_n(i,j)=(y0-tand(A)*(x0)-y_n(i-1,j)+tand(B)*(x_n(i-1,j)))/(tand(B)-tand(A));
                y_n(i,j)=(tand(B)*(x_n(i,j)-x_n(i-1,j)))+y_n(i-1,j);
            end
        elseif j>1
            if i==1 % Centerline
                A=((theta_n(i+1,j-1)+theta_n(i,j))/2)-((mu_n(i+1,j-1)+mu_n(i,j))/2);
                B=0;
                x_n(i,j)=(y_n(i+1,j-1)-tand(A)*(x_n(i+1,j-1)))/(tand(B)-tand(A));
                y_n(i,j)=0;
            else %(i>1) Rest of the points
                A=((theta_n(i+1,j-1)+theta_n(i,j))/2)-((mu_n(i+1,j-1)+mu_n(i,j))/2);
                B=((theta_n(i-1,j)+theta_n(i,j))/2)+((mu_n(i-1,j)+mu_n(i,j))/2);
                x_n(i,j)=(y_n(i+1,j-1)-tand(A)*(x_n(i+1,j-1))...
                    -y_n(i-1,j)+tand(B)*(x_n(i-1,j)))/(tand(B)-tand(A));
                y_n(i,j)=(tand(B)*(x_n(i,j)-x_n(i-1,j)))+y_n(i-1,j);
            end
        end
    end
end
% Initializing array for coordinates of the wall nodes
x_wall=zeros(1,n+1);
y_wall=zeros(1,n+1);
theta_wall=zeros(1,n+1);

% Point 'A' at the expansion point
x_wall(1,1)=x0;
y_wall(1,1)=y0;
theta_wall(1,1)=theta_max;

% Computation of the coordinates of rest of the wall nodes
for j=2:n+1
    if j==2
        A=(theta_wall(1,1));
        B=((theta_n(n-j+2,j-1)))+((mu_n(n-j+2,j-1)));
        x_wall(1,j)=(y_wall(1,j-1)-tand(A)*(x_wall(1,j-1))...
            -y_n(n-j+2,j-1)+tand(B)*(x_n(n-j+2,j-1)))/(tand(B)-tand(A));
        y_wall(1,j)=tand(theta_wall(1,1))*(x_wall(1,j)-x_wall(1,j-1))+y_wall(1,j-1);
    elseif j>2
        A=(theta_n(n-j+3,j-2)+(theta_n(n-j+2,j-1)))/2;
        B=((theta_n(n-j+2,j-1))+(mu_n(n-j+2,j-1)));
        x_wall(1,j)=(y_wall(1,j-1)-tand(A)*(x_wall(1,j-1))...
            -y_n(n-j+2,j-1)+tand(B)*(x_n(n-j+2,j-1)))/(tand(B)-tand(A));
        y_wall(1,j)=tand(A)*(x_wall(1,j)-x_wall(1,j-1))+y_wall(1,j-1);
    end
end

%% Plotting of the nozzle contour and characteristic lines
figure(1)
plot(x_wall,y_wall,'b-','LineWidth',2)
yline(0,'k--');
xline(0,'k--');
grid on
hold on
% Plotting the C- Characteristics from 'a' to points on First Characteristic line 
for j=1:n 
 plot([x0 x_n(j,1)],[y0 y_n(j,1)],'b-o','MarkerFaceColor','k','MarkerSize',3)
end
% Plotting the C+ Characteristics  
for j=1:n-1
     plot(x_n(1:n-j+1,j),y_n(1:n-j+1,j),'b-o','MarkerFaceColor','k','MarkerSize',3)
end
% Plotting and connecting the Inner nodes  
 for j=1:n
    for i=2:n-j+1
        plot([x_n(i,j) x_n(i-1,j+1)],[y_n(i,j) y_n(i-1,j+1)],'b-o','MarkerFaceColor','k','MarkerSize',3)
    end
 end
% Joining the C+ Characteristics with Wall nodes  
for i=2:n+1
        plot([x_wall(1,i) x_n(n-i+2,i-1)],[y_wall(1,i) y_n(n-i+2,i-1)],'b-o','MarkerFaceColor','k','MarkerSize',3.5)
end

title('MLN Contour with Method of Characteristics','fontsize',14,'interpreter','latex')
xlabel('$x$','fontsize',14,'interpreter','latex')
ylabel('$y$','fontsize',14,'interpreter','latex')
xlim([x0-0.5 max(x_wall)])
ylim([y0-1.0 max(y_wall)])
axis equal

%% Writing the Data points of the MLN Contour
x_contour_data= [x_wall max(x_wall) 0 x_wall(1,1)];
y_contour_data= [y_wall 0 0 y_wall(1,1)];
num=length(x_contour_data);
z_contour_data=zeros(1,num);
Group_num=ones(1,num).*2;
a=1:num;
data=[x_wall;y_wall;zeros(1,length(x_wall))];
file=fopen('MLN_Contour_Soliworks_Data.txt','w');
fprintf(file,'%d %d %d\n',data);
fclose(file);

% Plot with Nozzle Datapoints
figure(2)
plot(x_contour_data,y_contour_data,'b-','LineWidth',2)
yline(0,'k--');
xline(0,'k--');
grid on
title('MLN Contour with Method of Characteristics','fontsize',14,'interpreter','latex')
xlabel('$x$','fontsize',14,'interpreter','latex')
ylabel('$y$','fontsize',14,'interpreter','latex')
xlim([x0-0.5 max(x_wall)])
ylim([y0-1.0 max(y_wall)])
axis equal

