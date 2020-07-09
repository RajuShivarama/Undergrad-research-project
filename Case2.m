%---------------------Input variables---------------------%
c=0.5;                           %----Radial clearance in mm
L=70;                            %----Length of bearing in mm
R=35;                            %----Radius of journal in mm
e=0.3;                           %----Eccentricity of shaft in mm
s=8;                             %----Bump pitch in mm
thick=0.25;                          %----Bump foil thickness in mm
pa=0.1;                          %----Atmospheric pressure in N/mm^2
v=0.3;                           %----Poissions Ratio 
n1=0.63*10^-6;                   %Viscosity of lubricant(MR fluid) in Ns/mm^2
w=2100/60;                        %Angular speed of rotation of shaft in rps 
Eb=210000;                       %Modulas of Elasticiy of bump foil in Mpa
alpha=(2*pa*s*(1-v^2))/(c*Eb*thick^3);%Compliance number of bump foil bearing 
Er=e/c;                          %--Ecentricity Ratio
B1=(6*w*n1/pa)*(R/c)^2;           %----Compressibility or Bearing number
 
%---------------------Matlab code by Central Finite differnce method-----------------------%
n=50;                            
m=50;                            
dt=(2*pi)/n;                     
dr=c/m;                        
iter=1000000;                  
l=linspace(0,70,50);      
r=linspace(0,c,50);
t=linspace(0,2*pi/50,50);
theta=2*pi/36;
r1=c;
 
for i=1:n
    for j=1:m
        p(i,j)=1.0;             
    end
end
sum(1)=1.0;
for k=1:iter
    sumij=0.0;
    
    for i=2:n-1  
        h=1+Er*cos(t)+(-2+alpha/c*(1.5*pa*p(i,j)-1));      %---h(i,j)
        h1=1+Er*cos(t)+alpha/c*(pa*p(i,j)-1)-0.5*dr;       %---h(i-0.5,j) 
        h2=1+Er*cos(t)+alpha/c*(pa*p(i,j)-1)+0.5*dr;       %---h(i+0.5,j) 
        h3=1+Er*cos(t)+alpha/c*(pa*p(i,j)-1)-dr;           %---h(i-1,j) 
        h4=1+Er*cos(t)+alpha/c*(pa*p(i,j)-1)+dr;           %---h(i+1,j) 
        h5=1+Er*cos(t)+alpha/c*(pa*p(i,j)-1)-0.5*dt;       %---h(i,j-0.5)
        h6=1+Er*cos(t)+alpha/c*(pa*p(i,j)-1)+0.5*dt;       %---h(i,j+0.5)
        h7=1+Er*cos(t)+alpha/c*(pa*p(i,j)-1)-dt;           %---h(i,j-1) 
        h8=1+Er*cos(t)+alpha/c*(pa*p(i,j)-1)+dt;           %---h(i,j+1) 
       
        cubhm=h1.*h1.*h1;
        cubhp=h2.*h2.*h2;
        cubhm1=h5.*h5.*h5;
        cubhp2=h6.*h6.*h6;
        
        a1=theta^2*dr^2;
        a2=r1*theta*dt*dr;
        a3=r1^2*dt^2*100;
        a4=1/(theta*dr)+1/(r1*dt);
        a5=1/(theta*dr)-1/(r1*dt);
        const=(cubhp/dr+cubhp2/dt)*a4+(cubhm/dr+cubhm1/dt)*a5;
        A= (h2/a1+h6/a2)/const;
        B= (h2/a2+h6/a3)/const ;
        C= (h1/a1+h5/a2)/const;
        D= (h1/a2+h5/a3)/const;
        E=  ((h3-h4)/(2*theta*dr)+(h7-h8)/(2*r1*dt))/const;
        
         for j=2:m-1;
             
          p(i,j)=A*p(i+1,j)+B*p(i,j+1)+A*p(i+1,j)+C*p(i-1,j)+D*p(i,j-1)+E;                                                
             sumij=sumij+p(i,j);
            
         end
    end
    sum(k+1)=sumij;
    percentage=abs(sum(k+1)-sum(k))/abs(sum(k+1));                
    if percentage<0.00001
        break
    end
end
y=k;
P=max(p);
%----------------Differnt Plots of pressure profile,Film thickness and Load carrying capacity---------------------% 
figure(1);                                              
surf(p)                     
figure(2);                                                
[X,Z]=meshgrid(t,r);
mesh(t,r,p);
title('Pressure profile of journal Foil bearing ');
xlabel('along one bump(2*pi/50)');
zlabel(' pressure variation');
ylabel('along the radial clearnce')
 
figure(3);                                                  
plot(t,P,'-k')  
title('Pressure distribution in 2D along circumferential');
xlabel('Theta in radian along one bump');
ylabel('pressure variation in N/mm2');
 
figure(6);                                                  
plot(t,h,'-k');
title('The variation of film thickness in circuferential direction')
xlabel('Theta in radian along one bump');
ylabel('film thickness in mm)');   
