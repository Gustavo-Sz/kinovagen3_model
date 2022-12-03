%% questao 5 e 6
%% jacobinano com origem do sistema 7 deslocado para coincidir 
% com do sistema 6
conf_1 = [0,pi/2,0,-pi/2,0,pi/4,0]
jacob0_1p = gen3_punho.jacob0(conf_1)
j_1p= round(j_p(0,pi/2,0,-pi/2,0),4)

conf_2 = [0,pi/4,0,pi/2,0,-pi/4,0]
jacob0_2p = gen3_punho.jacob0(conf_2)
j_2p = round(j_p(0,pi/4,0,pi/2,0),4)

conf_3 = [0,pi/2,0,-pi/4,0,pi/2,0]
jacob0_3p = gen3_punho.jacob0(conf_3)
j_3p= round(j_p(0,pi/2,0,-pi/4,0),4)

%% questao 7
%% novo jacobiano considerando sistema de coordenadas 7 no efetuador, 
% como originalmente estava.

% validacao
conf_1 = [0,pi/2,0,-pi/2,0,pi/4,0]
j1 = gen3.jacob0(conf_1)
j1_m = round(j_t(0,pi/2,0,-pi/2,0,pi/4,0,-167.4/1000),4)

conf_2 = [0,pi/4,0,pi/2,0,-pi/4,0]
j2 = gen3.jacob0(conf_2)
j2_m = round(j_t(0,pi/4,0,pi/2,0,-pi/4,0,-167.4/1000),4)

conf_3 = [0,pi/2,0,-pi/4,0,pi/2,0]
j3 = gen3.jacob0(conf_3)
j3_m = round(j_t(0,pi/2,0,-pi/4,0,pi/2,0,-167.4/1000),4)


% calculo dos torques
jac = j_t(0,15*(pi/180),pi,230*(pi/180),0,55*(pi/180),pi/2,-167.4/1000);
torq_analitico = transpose(round(jac,5))*[0;0;-25;0;0;0]
torq_model = gen3.jacob0([0,15*(pi/180),pi,230*(pi/180),0,55*(pi/180),pi/2])'*[0;0;-25;0;0;0]

function j = j_p(t1,t2,t3,t4,t5)

function t = t_ab(ang,alp,a,d)
t = [cos(ang) -sin(ang)*round(cos(alp)) sin(ang)*round(sin(alp)) a*cos(ang);
    sin(ang) cos(ang)*round(cos(alp)) -cos(ang)*round(sin(alp)) a*sin(ang);
    0 round(sin(alp)) round(cos(alp)) round(d,3); 0 0 0 1];
end
t01 = t_ab(t1,pi/2,0,-284.8/1000);
t12 = t_ab(t2,-pi/2,0,-11.8/1000);
t23 = t_ab(t3,pi/2,0,-420.8/1000);
t34 = t_ab(t4,-pi/2,0,-12.8/1000);
t45 = t_ab(t5,pi/2,0,-314.3/1000);

r01 = t01(1:3,1:3);
r12 = t12(1:3,1:3);
r23 = t23(1:3,1:3);
r34 = t34(1:3,1:3);
r45 = t45(1:3,1:3);

p12 = t12(1:3,4);
p23 = t23(1:3,4);
p34 = t34(1:3,4);
p45 = t45(1:3,4);

r02 = r01*r12;
r03 = r02*r23;
r04 = r03*r34;
r05 = r04*r45;

rb0 = roty(180);
rb1 = rb0*r01;
rb2 = rb0*r02;
rb3 = rb0*r03;
rb4 = rb0*r04;
rb5 = rb0*r05;

p55_b = rb5*[0;0;0];
p45_b = rb4*p45 + p55_b;
p35_b = rb3*p34 + p45_b;
p25_b = rb2*p23 + p35_b;
p15_b = rb1*p12 + p25_b;

h1_b = rb0*[0;0;1];
h2_b = rb1(1:3,3);
h3_b = rb2(1:3,3);
h4_b = rb3(1:3,3);
h5_b = rb5(1:3,3);

j1_13 = skew(h1_b)*p15_b;
j2_13 = skew(h2_b)*p25_b;
j3_13 = skew(h3_b)*p35_b;
j4_13 = skew(h4_b)*p45_b;
j5_13 = skew(h5_b)*p55_b;

j = [j1_13 j2_13 j3_13 j4_13 j5_13];
end

function j = j_t(t1,t2,t3,t4,t5,t6,t7,d7)

function t = t_ab(ang,alp,a,d)
t = [cos(ang) -sin(ang)*round(cos(alp)) sin(ang)*round(sin(alp)) a*cos(ang);
    sin(ang) cos(ang)*round(cos(alp)) -cos(ang)*round(sin(alp)) a*sin(ang);
    0 round(sin(alp)) round(cos(alp)) round(d,3); 0 0 0 1];
end
t01 = t_ab(t1,pi/2,0,-284.8/1000);
t12 = t_ab(t2,-pi/2,0,-11.8/1000);
t23 = t_ab(t3,pi/2,0,-420.8/1000);
t34 = t_ab(t4,-pi/2,0,-12.8/1000);
t45 = t_ab(t5,pi/2,0,-314.3/1000);
t56 = t_ab(t6,-pi/2,0,0);
t67 = t_ab(t7,pi,0,d7);

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

r02 = r01*r12;
r03 = r01*r12*r23;
r04 = r01*r12*r23*r34;
r05 = r01*r12*r23*r34*r45;
r06 = r01*r12*r23*r34*r45*r56;

rb0 = roty(180);
p77_b =[0;0;0];
p67_b = rb0*r06*p67;
p57_b = rb0*r05*p56 + p67_b;
p47_b = rb0*r04*p45 + p57_b;
p37_b = rb0*r03*p34 + p47_b;
p27_b = rb0*r02*p23 + p37_b;
p17_b = rb0*r01*p12 + p27_b;

h1_b = rb0*[0;0;1];
h2_b = rb0*r01(1:3,3);
h3_b = rb0*r02(1:3,3);
h4_b = rb0*r03(1:3,3);
h5_b = rb0*r04(1:3,3);
h6_b = rb0*r05(1:3,3);
h7_b = rb0*r06(1:3,3);

j1_13 = skew(h1_b)*p17_b;
j2_13 = skew(h2_b)*p27_b;
j3_13 = skew(h3_b)*p37_b;
j4_13 = skew(h4_b)*p47_b;
j5_13 = skew(h5_b)*p57_b;
j6_13 = skew(h6_b)*p67_b;
j7_13 = skew(h7_b)*p77_b;

j1_36 = h1_b;
j2_36 = h2_b;
j3_36 = h3_b;
j4_36 = h4_b;
j5_36 = h5_b;
j6_36 = h6_b;
j7_36 = h7_b;
j = [j1_13 j2_13 j3_13 j4_13 j5_13 j6_13 j7_13; j1_36 j2_36 j3_36 j4_36 j5_36 j6_36 j7_36];

end 
