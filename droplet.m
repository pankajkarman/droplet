
%%
R = (0.75*0.001*.5)/(cos(pi/10));
rho=925;
k=2.219;
Cp=(2.1*1000);
alpha=10;
Ta=10;
Tc=-10;
Tc_l=.01;
delta_t=.01;
rho_l=998.2;
k_l = 0.6;
Cp_l = (4.18*1000);
Qv=333000;
M=((rho_l)*(2*pi*R^3)/3);
delta_r=(R/48);
delta_q=(pi/20);
r1=0;
q1=0;
S=(alpha*delta_r/k);
S_l=(alpha*delta_r/k_l);
U=(k*delta_t)/(rho*Cp);
U_l=(k_l*delta_t)/(rho_l*Cp_l);
Ts=283000;
%P=U_l/(delta_r)^2;
%L=cot(delta_q);
%%
%%
%initial condition
T(1:49,1,1:Ts)=Tc;
T(1:49,21,1:Ts)=Tc;
T(2:49,2:20,1)=Ta;

%%
for n=1:Ts
    for i=2:48
        for j=2:20
            %boundary condition
         
               
            
        
        
           
            
 T(i,j,n+1)=T(i,j,n)+(U_l*(1/((r1+((i-1)*delta_r))*(delta_q)^2))*(T(i,j+1,n)+T(i,j-1,n)-(2*T(i,j,n))));
 
    if T(i,j,n+1)==0
        
         T(i,j,n+2)=T(i,j,n+1)+(U*(1/(r1+((i-1)*delta_r))*(delta_q)^2)*(T(i,j+1,n+1)+T(i,j-1,n+1)-(2*T(i,j,n+1))))+((Qv*delta_t)/(rho*Cp));
    
 
  elseif T(i,j,n+1)<0
        
  T(i,j,n+2)=T(i,j,n+1)+(U*(1/((r1+((i-1)*delta_r))*(delta_q)^2))*(T(i,j+1,n+1)+T(i,j-1,n+1)-(2*T(i,j,n+1))));
 
    end
   
    
    if (T(i,10,n)>-(0.001)&& T(i,10,n)<.001)
        m1(n)=i*delta_r;
        m=m1(n);
        b(n)=sqrt(abs((R^2)-(m^2))) ;
        
        a=b(n);
        M_11(n)=rho*pi*(((R^2)*m)-(m^3/3));
        M_1=M_11(n);
        M_22(n)=M-M_1 ;
        M_2=M_22(n);
        P1(n)=sqrt(((6*M_2)/(2*pi*rho_l))^2+a^3);
        P=P1(n);
        A1(n)=(((6*M_2)/(2*pi*rho_l))+P);
        A=A1(n);
        B1(n)=(((6*M_2)/(2*pi*rho_l))-P);
        B=B1(n);
        g(n)=nthroot(A,3)+nthroot(B,3);
       
        h=g(n);
       
       
     
    end
        end
    end
end

%%



temp=T(:,:,30000);
%%
q=linspace(0,pi,21);
r=linspace(0,R,49);

[q,r]=meshgrid(q,r);
[x,y]=pol2cart(q,r);

xx=0.01;
yy=0.01;

reruns=1;                  % number of times movie is to play
fps=500;                     % frames per second

nframes = Ts;              % number of frames in the movie
Frames = moviein(nframes);
for i=1:2000:Ts
contourf(x,y,T(:,:,i))
colormap gray
colorbar('Eastoutside')
axis equal
title(i);
Frames(:,:) = getframe;
end
 movie(Frames,reruns,fps)
 %%
 h3=figure;
 plot(b,g,'+');
 title('b-g');
h4=figure;
 plot(b,m1,'+');
 title('b-m1');
 
 h6=figure;
 plot(b,M_11,'+');
 title('b-M_11');
 h7=figure;
 plot(b,M_22,'+');
 title('b-M_22');
 f8=figure;
 plot(b,A1,'+');
 title('b-A1');
 
 h9=figure;
 plot(b,B1,'+');
 title('b-B1');
 C=(6*M_22)/(2*pi*rho_l);
 h10=figure;
 plot(b,C,'+');
 title('b-c')
hold on
 plot(b,P1,'o');
 legend('C','P')
 hold off

        
 
 
 
 