%% QUESTAO 2 %%
fac= 1/1000;
d1 = fac*284.8;
d2 = fac*11.8;
d3 = fac*420.8;
d4 = fac*12.8;
d5 = fac*314.3;
d7 = fac*167.4;

%All link lengths and offsets are measured in m 
clear links
%            theta      d           a       alpha
links = [
	    Link([0        -d1           0       pi/2   0], 'standard')
		Link([0        -d2           0       -pi/2  0], 'standard')
		Link([0        -d3           0       pi/2   0], 'standard')
		Link([0        -d4           0       -pi/2  0], 'standard')
		Link([0        -d5           0       pi/2   0], 'standard')
		Link([0         0            0       -pi/2  0], 'standard')
        Link([0        -d7           0       pi     0], 'standard')

	];

gen3=SerialLink(links, 'name', 'Kinova Gen3');

% o sistema de coordenadas da base do braço e tal que a transformação para
% o sistema de coordenadas do primeiro elo é uma rotação de 180° em torno
% de x, ou seja o eixo Z do sistema da base está em sentido oposto ao eixo
% Z dos elos que estão na vertical, considerando a configuração inicial.
gen3.base=troty(180);

%%

%% QUESTAO 3 %% Validação de parâmetros de DH Standard
qz = [0 0 0 0 0 0 0];
qd = [0 41.4 0 90 0 -41.4 0]*pi/180;
qr = [0 0  0 90 180 0 0]*pi/180;
q1 = [0 90 0 -90 0 90 0]*pi/180;

"q1"
t1= t_b7(0, pi/2, 0, -pi/2, 0, pi/2, 0)
t1_m = gen3.fkine(q1)
"qr"
t2=t_b7(0, 0, 0, pi/2, pi, 0, 0)
t2_m = gen3.fkine(qr)
"qd"
t3=t_b7(0, 41.4*(pi/180), 0, pi/2, 0, -41.4*(pi/180), 0)
t3_m = gen3.fkine(qd)
"qz"
t4=t_b7(0, 0, 0, 0, 0, 0, 0)
t4_m = gen3.fkine(qz)
%% 

%% 	QUESTAO 4 %% Verificar cinemática inversa
r = [0 0 1; 1 0 0; 0 1 0];

pa = [0.58; 0; 0.43];
pb = [0.70;0;0.1]; 
pc = [0.78;0;0.10];
pd = [0.58;0;0.43];
pe = [0.58;-0.15;0.43];
pf = [0.78;-0.15;0.1];
pg = [0.7;-0.15;0.1];
ph = [0.58;0;0.43];
p=[pa pb pc pd pe pf pg ph];

"pa"
tr_inv(gen3,r,pa);
"pb"
tr_inv(gen3,r,pb);
"pc"
tr_inv(gen3,r,pc);
"pd"
tr_inv(gen3,r,pd);
"pe"
tr_inv(gen3,r,pe);
"pf"
tr_inv(gen3,r,pf);
"pg"
tr_inv(gen3,r,pg);
"ph"
tr_inv(gen3,r,ph);

inv_explore(gen3,r,p);
%%

function [t, q_inv, t_inv,err] = tr_inv(gen,r,p)
    t = rt2tr(r,p);
    q_inv = gen.ikine(t);
    t_inv = gen.fkine(q_inv);
    err= sqrt(sum((t_inv.t - p).^2));
end

function [err] = inv_explore(gen,r,p)
    tol = [1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4];
    for i = 1:size(tol,2)
        fprintf("tol: %i\n",tol(i))
        mean = 0;
        tic;
        for j = 1:size(p,2)
            t = rt2tr(r,p(:,j));
            q_inv = gen.ikine(t,'tol',tol(i));
            t_inv = gen.fkine(q_inv);
            %fprintf("p %i: %f,%f,%f\n",j, p(1,j),p(2,j),p(3,j));
            mean = mean + sqrt(sum((t_inv.t - p(:,j)).^2));
        end
        time = toc
        time = 0;
        mean = mean/size(p,2)
    end
    q0 = [[0 90 0 0 0 0 0]; [90 90 0 0 0 0 0 ]; [180 90 0 0 0 0 0];[180 65 0 120, 0 90 0] ]*(pi/180)
    for i = 1:size(q0,1)
        mean = 0;
        tic;
        fprintf("q0: q%i\n",i)
        for j = 1:size(p,2)
            t = rt2tr(r,p(:,j));
            q_inv = gen.ikine(t,q0(i,:),'tol',tol(i));
            t_inv = gen.fkine(q_inv);
            %fprintf("p %i: %f,%f,%f\n",j, p(1,j),p(2,j),p(3,j));
            mean = mean + sqrt(sum((t_inv.t - p(:,j)).^2));
        end
        time = toc
        time = 0;
        mean = mean/size(p,2)
    end

end


function tb7 = t_b7(t1,t2,t3,t4,t5,t6,t7) 

function t = t_ab(t,alf,a,d) 
    t = [cos(t) -sin(t)*round(cos(alf)) sin(t)*round(sin(alf)) a*cos(t);
        sin(t) cos(t)*round(cos(alf)) -cos(t)*round(sin(alf)) a*sin(t);
        0 round(sin(alf)) round(cos(alf)) round(d,3); 
        0 0 0 1];
end

t01 = t_ab(t1,  pi/2, 0, -284.8/1000);
t12 = t_ab(t2, -pi/2, 0, -11.8/1000);
t23 = t_ab(t3,  pi/2, 0, -420.8/1000);
t34 = t_ab(t4, -pi/2, 0, -12.8/1000);
t45 = t_ab(t5,  pi/2, 0, -314.3/1000);
t56 = t_ab(t6, -pi/2, 0, 0);
t67 = t_ab(t7,  pi,   0, -167.4/1000);

r01 = t01(1:3,1:3);
r12 = t12(1:3,1:3);
r23 = t23(1:3,1:3);
r34 = t34(1:3,1:3);
r45 = t45(1:3,1:3);
r56 = t56(1:3,1:3);
r67 = t67(1:3,1:3);

p01 = t01(1:3,4);
p12 = t12(1:3,4);
p23 = t23(1:3,4);
p34 = t34(1:3,4);
p45 = t45(1:3,4);
p56 = t56(1:3,4);
p67 = t67(1:3,4);

rb0 = roty(180);
r02 = r01*r12;
r03 = r01*r12*r23;
r04 = r01*r12*r23*r34;
r05 = r01*r12*r23*r34*r45;
r06 = r01*r12*r23*r34*r45*r56;
r07 = r01*r12*r23*r34*r45*r56*r67;
rb7 = rb0*r07;

pb7 = rb0*p01 + rb0*r01*p12 + rb0*r02*p23 +rb0*r03*p34 + rb0*r04*p45 + rb0*r05*p56 + rb0*r06*p67;
tb7 = [rb7 pb7;0 0 0 1 ];

end
