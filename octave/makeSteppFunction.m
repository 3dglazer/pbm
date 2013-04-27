function [samples,values] =makeSteppFunction (x,y)
	dPoints=x;
	n=max(size(x));
	vals=y;
	steppedPoints=zeros(1,n*2);
	steppedVals=zeros(1,n*2);
	for i=[1:n-1]
		id1=((i-1)*2+1);
		id2=((i-1)*2+2);
		steppedPoints(id1)=dPoints(i);
		steppedPoints(id2)=dPoints(i+1);
		steppedVals(id1)=vals(i);
		steppedVals(id2)=vals(i);
	end
	steppedPoints(n*2-1)=dPoints(n);
	steppedPoints(n*2)=dPoints(n);
	steppedVals(n*2-1)=vals(n);
	steppedVals(n*2)=vals(n);
	samples=steppedPoints;
	values=steppedVals;
end
