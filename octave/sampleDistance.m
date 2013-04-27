function [dist] = sampleDistance (x,omega,sigmaTMax,funname)
 s=0;
 maxIt=1000;
 it=0;
 while(it<maxIt)
	s+=-(log(1.-rand())/sigmaTMax);
	%misto zavorky by melo byt sigmaT(x) pro nehomog prostredi
	movedX=x+s*omega;
	sigmaTX=feval (funname, movedX);
	probab=sigmaTX/sigmaTMax;
	if rand()<probab
		dist=s;
		return;
	end
	it+=1;
 end
 printf("\nMax iterations reached in woodcock!!!\n\n")
end