function [vals] = functionSum(evalPoints,funname)
sz=size(evalPoints);
sz=max(sz(1),sz(2));
vals=zeros(sz);
vsum=0;
for i=1:sz
	vsum+=feval(funname,evalPoints(i));
	vals(i)=vsum;
end
end
