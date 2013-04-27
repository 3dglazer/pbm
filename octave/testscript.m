close all
clear all

num=4;
funTransmition="myLinear";
funTransmition="myConstant";
maxSigma=0.5
[ds,valsWoodcock]=takeNWoodcockSamples (num,maxSigma,funTransmition);

%ploting the original func and corresponding transmitance
nEvals=300;
if (num>nEvals)
nEvals=num*2;
end

%ve vals je opticka tloustka nascitana pres vzdalenost musim ji predelat na transmitance, kvuli zpravne opticke tlousce musi sedet zacatky puvodni funkce a samplu
evalPoints=linspace(min(ds),max(ds),nEvals);
vals=opticalThickness(evalPoints,funTransmition);
vals=exp(-1*vals); %predela to na transmitanci
figure
axis equal

%ploting transmitance
plot(evalPoints,vals,'r');

%ploting the woodcock samples
hold on
[ds,vals]=makeSteppFunction(ds,valsWoodcock);
plot(ds,vals,'b');
hold on
plot(ds,vals,'oy');
