fac= 1/1000;
d1 = fac*284.8;
d2 = fac*11.8;
d3 = fac*420.8;
d4 = fac*12.8;
d5 = fac*314.3;

%All link lengths and offsets are measured in m
clear links
%            theta    d           a       alpha
links = [
	    Link([0        -d1         0       pi/2  0], 'standard')
		Link([0        -d2         0       -pi/2 0], 'standard')
		Link([0        -d3         0       pi/2  0], 'standard')
		Link([0        -d4         0       -pi/2 0], 'standard')
		Link([0        -d5         0       pi/2  0], 'standard')
		Link([0        0           0       -pi/2 0], 'standard')
        Link([0        0           0       pi    0], 'standard')

	];

gen3_punho=SerialLink(links, 'name', 'Kinova Gen3');
gen3_punho.base=troty(180);