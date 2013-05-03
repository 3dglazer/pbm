close all
clear all

num=16;
funTransmition="myLinear";
funTransmition="myConstant";
maxSigma=0.5
%ploting the original func and corresponding transmitance
nEvals=600;
if (num>nEvals)
nEvals=num*2;
end
%will progressively refine the estimate
progressiveTransmittance=zeros(1,nEvals);
densEval=linspace(0,600,nEvals);
%how many itterations
niter=0;
for i=[1:16]
niter+=1;
[ds,valsWoodcock]=takeNWoodcockSamples (num,maxSigma,funTransmition);

prgTr=interp1(ds, valsWoodcock, densEval,'nearest');
fineMask=finite(prgTr); %filters out NANs
prgTr(isnan(prgTr)) = 0; %sets NANs to zero
progressiveTransmittance(fineMask)+=prgTr(fineMask);
progressiveTransmittance(fineMask)*=0.5;
%ve vals je opticka tloustka nascitana pres vzdalenost musim ji predelat na transmitance, kvuli zpravne opticke tlousce musi sedet zacatky puvodni funkce a samplu
evalPoints=linspace(min(ds),max(ds),nEvals);
valsTr=opticalThickness(evalPoints,funTransmition);
valsTr=exp(-1*valsTr); %predela to na transmitanci

figure
axis equal

%ploting transmitance
plot(evalPoints,valsTr,'r');

%ploting the woodcock samples
hold on
[ds,vals]=makeSteppFunction(ds,valsWoodcock);
plot(ds,vals,'b');
hold on
plot(ds,vals,'oy');
end

figure
%ploting transmitance
plot(evalPoints,valsTr,'r');
hold on
plot(densEval,progressiveTransmittance,'g');