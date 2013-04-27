function [dist] = sampleDistance (x,omega,sigmaTMax)
 s=0;
 maxIt=10;
 it=0;
 while(it<maxIt)
	s+=-(log(1.-rand())/sigmaTMax);
	rn=rand();
	%misto zavorky by melo byt sigmaT(x) pro nehomog prostredi
	probab=(x+s*omega)*sigmaTMax;
	if rn<probab
		dist=s;
		break;
	end
	it+=1;
 end
end