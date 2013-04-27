function [samples,values] =takeNWoodcockSamples (n,maxSigma,funcname)
	dPoints=zeros(1,n+1);
	funTransmition=funcname;
	printf("there is %d samples to be generated\n",n)
	for i=[1:n],
		d=sampleDistance(0,1.,maxSigma,funTransmition);
		dPoints(i)=d;
	end
	dPoints=sort(dPoints);
	values=opticalThickness(dPoints,funTransmition);
	values=exp(-1*values); %predela to na transmitanci
	samples=dPoints;
end
