% Main MC

clc; clear all; close all;
% arbitrary parameters

 x=[1.85255251645473,0.484459074644468,2.45141037308471,1.24155168491416,0.292605862266227,2.44950898784127,4.77213132545400,0.158172824056465,2.16182272907231,0.248765927014658,0.0122546337943833,0.0357260727168561,0.272970092608579,1.92804613252871,0.317398440514774,0.108703248928418,0.0911836618960659,1.15999178783990,0.0130815489973584,-0.898507677965659,-0.0793197223439565,0.00391616510303808,0.0646978829533495,0.0137580598397501,0.5];


 NTime=250; NSim=100000;T=1;T0=1;F0=100;t=0;
 
 K=80:5:120; %Strikes

  
 price_rough = option_rough_forward_mc(x,NTime,NSim,t, T0, T,K,F0) % Rough forward model (if H=0.5, then these prices should be the same as for the classical model).

 price_classical = Option_classical_mc(x,t,T0,T,F0,K, NSim, NTime) % Classical forward model 