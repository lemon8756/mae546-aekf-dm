function out = MakeNewRippleMask(filename,N)

xs = (-N+0.5:N-0.5);
ys = xs';

TMP=load(filename);
L=length(TMP);
Mask=zeros(2*L,2);

   Mask(1:80,1)=TMP(1:80,1);
   Mask(1:80,2)=TMP(1:80,2);
   Mask(81:160,1)=TMP(1:80,1);
   Mask(81:160,2)=TMP(1:80,3);
   
   Mask(161:240,1)=TMP(81:160,1);
   Mask(161:240,2)=TMP(81:160,2);
   Mask(241:320,1)=TMP(81:160,1);
   Mask(241:320,2)=TMP(81:160,3);
   
   Mask(321:400,1)=TMP(161:240,1);
   Mask(321:400,2)=TMP(161:240,2);
   Mask(401:480,1)=TMP(161:240,1);
   Mask(401:480,2)=TMP(161:240,3);
   
   Mask(481:560,1)=TMP(241:320,1);
   Mask(481:560,2)=TMP(241:320,2);
   Mask(561:640,1)=TMP(241:320,1);
   Mask(561:640,2)=TMP(241:320,3);

out = zeros(2*N,2*N);
sgn = -1;
for k=1:8
    A = Mask(1+(k-1)*80:k*80,:);
    [m, n] = size(A);
    A=[(A(m:-1:1,:)*[-1 0; 0 1])' A']';

    yofx = spline(2*N*A(:,1),2*N*A(:,2),xs);
    for j=1:2*N
	m1 = (ys < yofx(j)-0.5) .* (ys >= -yofx(j)+0.5);
	m2 = (ys < yofx(j)+0.5) .* (ys >= -yofx(j)-0.5);
	m21 = m2-m1;
        out(:,j) = out(:,j) + sgn*(m21*(yofx(j)-floor(yofx(j))) + m1);
    end
    sgn = -sgn;
end
