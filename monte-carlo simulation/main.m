%Made with MATLAB2022a
%% parameters
% generate (inside the s_measure program - gen) - if you want to test the programm set it to true
% parameters required to work with Monte Carlo simulation:
% run - number of Monte Carlo evolution runs, also used as the seed for the
% first run
% r - recombination probability per genome
% M - crossover number
% L - number of loci
% N - population size
% run2 - the seed for generating the random distribution of s
% NUbs2 - N*Ub*s^2 (Ub - probability of the benefitial mutation)
% tf - end time of evolution
% s0 - the width of the uniform distribution of selection coefficient
% f0 - initial value of f
% ac*s0, bc*s0 - borders of uniform s distribution
% mu - mutation probability per site
% tsec1, tsec2, tsec3 - time of the first, second and third sequence sample
% parameters required to work both with real and generated data
% C - initial value of C in equation 1
% appr - method of curve approximation:
% 'poly' - by basic polynomials
% 'spline' - by cubic splines
% 'pchip' - by Piecewise Cubi Hermite Polynomial (PCHIP)
% 'test' - without any approximation and C finding (only for curves printing)
% apprR - additional parameter for 'poly' approximation - rate of
% polynomial

global s0 N tf f0 run2 ac bc tsec1 tsec2 tsec3 l rndmf0 var_t_f0 delay_freq delay_time population_num


generate=true;
rndmf0=0; %1 - включить случайное распределение f0 по сайтам, 0 - отключить
var_t_f0=0; %1 - включить отложенные последовательности, 0 - отключить
delay_freq=0.1; %какая часть популяции будет позже вступать в эволюцию
delay_time=6; % время после которого встраиваются новые геномы
population_num = 1; %количество взаимодействующих популяций
%run=5;
%r=0;
%L=40;
%M=2;
N=1000;
run2=21;
%mu=0.07/L;
%NUbs2=0.05;
tf=30;
    s0=0.1;
    l=0.15;
    f0=0.5;
    ac=-1;
    bc=1;
    %mu=NUbs2/N/s0^2/L;
    tsec1=5;
    tsec2=15;
    tsec3=30;
    tsecs=[tsec1,tsec2,tsec3];
    C=0;
    appr='test'; 
    apprR=4;

if generate
    %ri = 0;
    %M=0;
for M=[1,3,10]
    for ri=[5,50,100]
       for L=[40,100,400]
           for run=[5,10,50]
            for i =1:run
                mu=0.07/L;
                r=ri/100;
                [dat{i},s,~,~,~,~,~]=recomb_2022_test(r,s0,ac,bc,M,L,N,tf,f0,l,i,run2,mu); 
                Data{i}=dat{i};
                  
            end
            name= sprintf("r%dM%drun%dL%d.mat",ri,M,run,L);
            save(name,"Data",'s','tsecs');
            clear Data s mu dat
           end
       end
    end
end

%gen1{i}=dat{i}{1};
%gen2{i}=dat{i}{2};
%gen3{i}=dat{i}{3};
%S_{i}=s;


order=1:L;
%sdis = s_measure(gen1.',gen2.',gen3.',order,C,appr,apprR,generate);

else
    load('BinData.mat','data');
genome1 = data{1,1};
genome2 = data{1,2};
genome3 = data{1,3};
order = data{1,4};
   sdis = s_measure(genome1,genome2,genome3,order,C,appr,apprR,generate);
end

