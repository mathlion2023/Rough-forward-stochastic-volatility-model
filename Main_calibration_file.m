% This file is to run both the futures and options pricing routines
% Mesias Alfeus 2023
% clearing the workspace
clc; close all; clear all;
% Loading Data
load('CrudeOil_all_v2.mat')

%% Initial parameters Calibrat

A = []; b = []; Aeq = []; beq = [];  
 nd=1;be=5000;en=6000; % between 2006 to 2010 daily data (options and futures on Crude Oil)
Nn=40;

mod=2; % 1 for rough else classical
% Options for fmincon
  opts = optimoptions(@fmincon,'Display','iter','TolFun',1e-10,'TolCon',1e-10,'MaxFunctionEvaluations',3500,'UseParallel','Always','Algorithm','sqp');

if mod==1
    x0=[1.85255251645473,0.484459074644468,2.45141037308471,1.24155168491416,0.292605862266227,2.44950898784127,4.77213132545400,0.158172824056465,2.16182272907231,0.248765927014658,0.0122546337943833,0.0357260727168561,0.272970092608579,1.92804613252871,0.317398440514774,0.108703248928418,0.0911836618960659,1.15999178783990,0.0130815489973584,-0.898507677965659,-0.0793197223439565,0.00391616510303808,0.0646978829533495,0.0137580598397501,1];
 % Upper & lower bound
 lb = [zeros(18,1)+0.001; -0.95*ones(3,1); zeros(3,1)+0.001;0.5];
 ub = [5*ones(18,1); 0.95*ones(3,1); 3*ones(3,1);1.5];
else
 lb = [zeros(18,1)+0.001; -0.95*ones(3,1); zeros(3,1)+0.001];
 x0=[1.85255251645473,0.484459074644468,2.45141037308471,1.24155168491416,0.292605862266227,2.44950898784127,4.77213132545400,0.158172824056465,2.16182272907231,0.248765927014658,0.0122546337943833,0.0357260727168561,0.272970092608579,1.92804613252871,0.317398440514774,0.108703248928418,0.0911836618960659,1.15999178783990,0.0130815489973584,-0.898507677965659,-0.0793197223439565,0.00391616510303808,0.0646978829533495,0.0137580598397501];
 ub = [5*ones(18,1); 0.95*ones(3,1); 3*ones(3,1)];
end
    j=1;
 tic;
 for id=be:nd:en
     if j==1
         x0=x0;
     else
         x0=xOut(j-1,:);
     end

          [xOut(j,:),fval(j)] = fmincon(@(x)obj_fun(x,id,TTM_all,futures_all, options_all,mod,Nn),x0,A,b,Aeq,beq,lb,ub,[],opts);

     j=j+1
 end
 runtime=toc