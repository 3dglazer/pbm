function [vals] = opticalThickness(evalPoints,funname)
sz=size(evalPoints);
sz=max(sz(1),sz(2));
vals=zeros(sz);
vsum=0;
%first optical thickness is always zero
vals(1)=0;
for i=2:sz
	vsum+=(evalPoints(i)-evalPoints(i-1))*feval(funname,evalPoints(i));
	vals(i)=vsum;
end
end
