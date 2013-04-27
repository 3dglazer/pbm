clear all
close all
n=100;
m_pdf=[0, 0, 2,10,1, 0,20,2,2,2,5, 0]
ln=size(m_pdf,2)
smpls=linspace(0,ln-1,n);
smpls_val=interp1 (linspace(0,ln-1,ln), m_pdf, smpls);
figure
plot([0:ln-1],m_pdf,'b');
hold on
plot(smpls,smpls_val,'r');
plot(smpls,smpls_val,'kx');

totalSum=0;
for i=[1:n]
	totalSum+=smpls_val(i);
	smpls_val(i)=totalSum;
end
totalSum=1/totalSum;
smpls_val=smpls_val*totalSum;%normalizing
smpls=smpls/(ln-1);
figure

plot(smpls,smpls_val,'b');

figure
nrnd=20;
maxVal=max(m_pdf);
for i=[1:10]
	clf ();
	randomPoints=unifrnd (0,1,[1,nrnd]);
	samples=interp1(smpls_val,smpls,randomPoints);
	hold on
	for sgline=samples*(ln-1)
		plot([sgline,sgline],[0,maxVal],'r');
	end
	hold on
	plot([0:ln-1],m_pdf,'b');
	sleep(1);
end