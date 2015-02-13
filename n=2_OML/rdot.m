% This function takes an n=2 Lindblad superoperator and returns fM(r) and fm(r), the max 
% and min values of dr/dt for each r. The interval r in [0,1] is sampled at (samples+1) 
% points. The function also returns the intercept, rcrit, of fM(r), and two control
% vectors, nnM and nnm, 3d vectors with norm one, that achieve the fM and fm.
% The Lindblad superoperator is specified in terms of the eigenvalues of the real
% part of the GKS matrix (a1, a2, a3) and the imaginary parts of the
% off-diagonal components (b1, b2, b3). fM and fm are plotted if shouldplot=1.

function[rcrit, rdot, nnM, nnm] = rdot(a1,a2,a3,b1,b2,b3,samples, shouldplot)

% Safer to make very very small bj's exactly zero

if abs(b1)<5e-12
    b1=0;
end
if abs(b2)<5e-12
    b2=0;
end
if abs(b3)<5e-12
    b3=0;
end

% Initialize

bnorm=sqrt(b1^2+b2^2+b3^2);
rdot=zeros(samples+1,2); rcrit=0;
nnM=zeros(samples+1,3); nnm=zeros(samples+1,3);
rdot(1,1)=bnorm; rdot(1,2)=-bnorm;
if (bnorm > 0)
	nnM(1,1)=b1/bnorm; nnM(1,2)=b2/bnorm; nnM(1,3)=b3/bnorm;
	nnm(1,1)=-b1/bnorm; nnm(1,2)=-b2/bnorm; nnm(1,3)=-b3/bnorm;
end

% Ensure input data is acceptable

if a2>a1 || a3>a2 || a3<0 || 4*a1*a2*a3 < a1*b1^2 + a2*b2^2 + a3*b3^2 - 1E-8 || ...
		4*a1*a2 < b3^2 - 1E-8
    'bad input'
    [a1,a2,a3,b1,b2,b3]
    return
end

% For bnorm = 0, the LM equations can be solved exactly and no polynomials are
% involved. fM is always nonpositive, so r=0.

if bnorm == 0
    rdot(:,1)=-(a2+a3)*((0:samples)+1)/samples;
    rdot(:,2)=-(a1+a2)*((0:samples)+1)/samples;
    nnM(:,1)=ones(samples+1,1);    nnm(:,3)=ones(samples+1,1);
    
% For exactly one non-zero b component, the LM equations can be solved. There are
% two standard solutions, but also one or two others that may come into play.

elseif b1 == 0 && b2 == 0
	nnM(:,3)=sign(b3)*ones(samples+1,1);    nnm(:,3)=-sign(b3)*ones(samples+1,1);
	rdot(:,1)=bnorm-(a1+a2)*((0:samples)+1)/samples;
	rdot(:,2)=-bnorm-(a1+a2)*((0:samples)+1)/samples;
	rcrit=bnorm/(a1+a2);
    if a1 > a3 && abs(b3) <= 2*(a1-a3)
       for j= ceil(abs(b3)/2/(a1-a3)*samples):samples
           rdot(j+1,1)=bnorm^2/4/j*samples/(a1-a3)-(a2+a3)*j/samples;
           nnM(j+1,3)=abs(b3)/2/j*samples/(a1-a3);
           nnM(j+1,1)=sqrt(1-(nnM(j+1,3))^2);
           if rcrit > abs(b3)/2/(a1-a3)
		if rdot(j+1,1) <= 1E-14 && rdot(j,1) > 0  
		 rcrit = (rdot(j+1,1)*(j-1) - rdot(j,1)*j )/samples/(rdot(j+1,1) -rdot(j,1));
               end
           end
       end
    end
    
elseif b1 == 0 && b3 == 0
	nnM(:,2)=sign(b2)*ones(samples+1,1);    nnm(:,2)=-sign(b2)*ones(samples+1,1);
	rdot(:,1)=bnorm-(a1+a3)*((0:samples)+1)/samples;
	rdot(:,2)=-bnorm-(a1+a3)*((0:samples)+1)/samples;
	rcrit=bnorm/(a1+a3);
    if a2 > a3 && abs(b2) <= 2*(a2-a3)
       for j= ceil(abs(b2)/2/(a2-a3)*samples):samples
           rdot(j+1,2)=-bnorm^2/4/j*samples/(a2-a3)-(a1+a2)*j/samples;
           nnm(j+1,2)=-b3/2/j*samples/(a2-a3);
           nnm(j+1,3)=sqrt(1-(nnm(j+1,2))^2);
       end
    end
    if a1 > a2 && abs(b2) <= 2*(a1-a2)
       for j= ceil(abs(b2)/2/(a1-a2)*samples):samples
           rdot(j+1,1)=bnorm^2/4/j*samples/(a1-a2)-(a2+a3)*j/samples;
           nnM(j+1,2)=abs(b2)/2/j*samples/(a1-a2);
           nnM(j+1,1)=sqrt(1-(nnM(j+1,2))^2);
           if rcrit > abs(b2)/2/(a1-a2)
               if rdot(j+1,1) <= 1E-14 && rdot(j,1) > 0
                   rcrit = (rdot(j+1,1)*(j-1) - rdot(j,1)*j )/samples/(rdot(j+1,1) -rdot(j,1));
               end
           end
       end
    end
    
elseif b2 == 0 && b3 == 0
	nnM(:,1)=sign(b1)*ones(samples+1,1);    nnm(:,1)=-sign(b1)*ones(samples+1,1);
	rdot(:,1)=bnorm-(a2+a3)*((0:samples)+1)/samples;
	rdot(:,2)=-bnorm-(a2+a3)*((0:samples)+1)/samples;
	rcrit=bnorm/(a2+a3);
    if a1 > a3 && abs(b1) <= 2*(a1-a3)
       for j= ceil(abs(b1)/2/(a1-a3)*samples):samples
           rdot(j+1,2)=-bnorm^2/4/j*samples/(a1-a3)-(a1+a2)*j/samples;
           nnm(j+1,1)=-b1/2/j*samples/(a1-a3);
           nnm(j+1,3)=sqrt(1-(nnm(j+1,1))^2);
       end
    end
    
% If exactly two b components are nonzero, we have to solve a two or four
% degree polynomial. There is one other possible solution that must be
% checked afterwards.

elseif b3 == 0
    if a1 == a2
        nnM(:,1)=b1*ones(samples+1,1)/bnorm;    nnm(:,1)=-b1*ones(samples+1,1)/bnorm;
        nnM(:,2)=b2*ones(samples+1,1)/bnorm;    nnm(:,2)=-b2*ones(samples+1,1)/bnorm;
        rdot(:,1)=bnorm-(a1+a3)*((0:samples)+1)/samples;
        rdot(:,2)=-bnorm-(a1+a3)*((0:samples)+1)/samples;
        rcrit=bnorm/(a1+a3);
        if (a2 > a3 && bnorm/2/(a1-a3) <= 1)
            for j= ceil(bnorm/2/(a1-a3)*samples):samples
                rdot(j+1,2)=-bnorm^2/4/j*samples/(a1-a3)-(a1+a2)*j/samples;
                nnm(j+1,1)=-b1/2/j*samples/(a1-a3);
                nnm(j+1,2)=-b2/2/j*samples/(a1-a3);
                nnm(j+1,3)=sqrt(1-(nnm(j+1,1))^2-(nnm(j+1,2))^2);
            end
        end
    else
        for j=1:samples
            r=j/samples;
            c(1)=1;c(2)=-2*(a1+a2);c(3)=a1^2+4*a2*a1+a2^2-(b1^2+b2^2)/4/r^2;
            c(4)=-2*a1*a2*(a1+a2)+(a1*b2^2+a2*b1^2)/2/r^2;
            c(5)=a1^2*a2^2-(a1^2*b2^2+a2^2*b1^2)/4/r^2;
            rts=roots(c); realno=0;
            for k=1:4
                if abs(imag(rts(k)))<0.00001
                    realno=realno+1;
                end
            end
            rrts=zeros(realno,1);place=1;
            for k=1:4
                if abs(imag(rts(k)))<0.0001
                    rrts(place)=real(rts(k));
                    place=place+1;
                end
            end
            test=zeros(realno,1);
            for k=1:realno
                n=zeros(3,1);
                n(1)=b1/2/(rrts(k)-a1)/r;n(2)=b2/2/(rrts(k)-a2)/r;
                n=n/norm(n);
                test(k)=b1*n(1)+b2*n(2)+r*(-a3+a1*(n(1)^2-1)+a2*(n(2)^2-1));
            end
            [rdot(j+1,1),ind1]=max(test);
            [rdot(j+1,2),ind2]=min(test);
            nnM(j+1,1)=b1/2/(rrts(ind1)-a1)/r;nnM(j+1,2)=b2/2/(rrts(ind1)-a2)/r;
            nnm(j+1,1)=b1/(rrts(ind2)-a1)/r/2;nnm(j+1,2)=b2/2/(rrts(ind2)-a2)/r;
	    if a2 != a3
            	fstar = -b1^2/4/r/(a1-a3) + -b2^2/4/r/(a2-a3) - r*(a1+a2);
	    	if fstar < rdot(j+1,2) && r^2 >= b1^2/4/(a1-a3)^2 + b2^2/4/(a2-a3)^2
                	rdot(j+1,2) = fstar; 
                	nnm(j+1,1) = -b1/2/r/(a1-a3);
                	nnm(j+1,2) = -b2/2/r/(a2-a3);
                	nnm(j+1,3) = sqrt(1-(nnm(j+1,1))^2-(nnm(j+1,2))^2);
            	end
	    end
            nnM(j+1,1:3)=nnM(j+1,1:3)/norm(nnM(j+1,1:3));
	    nnm(j+1,1:3)=nnm(j+1,1:3)/norm(nnm(j+1,1:3));
            if rdot(j+1,1) <= 1E-14 && rdot(j,1) >= 0
               rcrit = (j - 1 + rdot(j,1)/(rdot(j,1)-rdot(j+1,1)))/samples; 
            end
        end
    end
elseif b2 == 0
    if a1 == a3
        nnM(:,1)=b1*ones(samples+1,1)/bnorm;    nnm(:,1)=-b1*ones(samples+1,1)/bnorm;
        nnM(:,3)=b3*ones(samples+1,1)/bnorm;    nnm(:,3)=-b3*ones(samples+1,1)/bnorm;
        rdot(:,1)=bnorm-(a1+a2)*((0:samples)+1)/samples;
        rdot(:,2)=-bnorm-(a1+a2)*((0:samples)+1)/samples;
        rcrit=bnorm/(a1+a2);
    else
        for j=1:samples
            r=j/samples;
            c(1)=1;c(2)=-2*(a1+a3);c(3)=a1^2+4*a3*a1+a3^2-(b1^2+b3^2)/4/r^2;
            c(4)=-2*a1*a3*(a1+a3)+(a1*b3^2+a3*b1^2)/2/r^2;
            c(5)=a1^2*a3^2-(a1^2*b3^2+a3^2*b1^2)/4/r^2;
            rts=roots(c); realno=0;
            for k=1:4
                if abs(imag(rts(k)))<0.00001
                    realno=realno+1;
                end
            end
            rrts=zeros(realno,1);place=1;
            for k=1:4
                if abs(imag(rts(k)))<0.0001
                    rrts(place)=real(rts(k));
                    place=place+1;
                end
            end
            test=zeros(realno,1);
            for k=1:realno
                n=zeros(3,1);
                n(1)=b1/2/(rrts(k)-a1)/r;n(3)=b3/2/(rrts(k)-a3)/r;
                n=n/norm(n);
                test(k)=b1*n(1)+b3*n(3)+r*(-a2+a1*(n(1)^2-1)+a3*(n(3)^2-1));
            end
            [rdot(j+1,1),ind1]=max(test);
            [rdot(j+1,2),ind2]=min(test);
            nnM(j+1,1)=b1/2/(rrts(ind1)-a2)/r;nnM(j+1,3)=b3/2/(rrts(ind1)-a3)/r;
            nnm(j+1,1)=b1/(rrts(ind2)-a2)/r/2;nnm(j+1,3)=b3/2/(rrts(ind2)-a3)/r;
            if a1>a2 && a2>a3
	    	fstar = b1^2/4/r/(a1-a2) - b3^2/4/r/(a2-a3) - r*(a1+a3);
	    	if fstar < rdot(j+1,2) && r^2 >= b1^2/4/(a1-a2)^2 + b3^2/4/(a2-a3)^2
                	rdot(j+1,2) = fstar; 
                	nnm(j+1,1) = b1/2/r/(a1-a2);
                	nnm(j+1,3) = -b3/2/r/(a2-a3);
                	nnm(j+1,2) = sqrt(1-(nnm(j+1,1))^2-(nnm(j+1,3))^2);
            	end
	    	if fstar > rdot(j+1,1) && r^2 >= b1^2/4/(a1-a2)^2 + b3^2/4/(a2-a3)^2
            	    rdot(j+1,1) = fstar; 
            	    nnM(j+1,1) = b1/2/r/(a1-a2);
            	    nnM(j+1,3) = -b3/2/r/(a2-a3);
            	    nnM(j+1,2) = sqrt(1-(nnm(j+1,1))^2-(nnm(j+1,3))^2);
            	end
	    end
            nnM(j+1,1:3)=nnM(j+1,1:3)/norm(nnM(j+1,1:3));
	    nnm(j+1,1:3)=nnm(j+1,1:3)/norm(nnm(j+1,1:3));
            if rdot(j+1,1) <= 1E-14 && rdot(j,1) >= 0
               rcrit = (j - 1 + rdot(j,1)/(rdot(j,1)-rdot(j+1,1)))/samples; 
            end
        end
    end
elseif b1 == 0
    if a2 == a3
        nnM(:,2)=b2*ones(samples+1,1)/bnorm;    
	nnm(:,2)=-b2*ones(samples+1,1)/bnorm;
        nnM(:,3)=b3*ones(samples+1,1)/bnorm;    
	nnm(:,3)=-b3*ones(samples+1,1)/bnorm;
        rdot(:,1)=bnorm-(a1+a2)*((0:samples)+1)/samples;
        rdot(:,2)=-bnorm-(a1+a2)*((0:samples)+1)/samples;
        rcrit=bnorm/(a1+a2);
        if (a1 > a2 && bnorm/2/(a1-a3) <= 1)
            for j= ceil(bnorm/2/(a1-a3)*samples):samples
                rdot(j+1,1)=bnorm^2/4/j*samples/(a1-a3)-(a2+a3)*j/samples;
                nnM(j+1,2)=b2/2/j*samples/(a1-a3);
                nnM(j+1,3)=b3/2/j*samples/(a1-a3);
                nnM(j+1,1)=sqrt(1-(nnM(j+1,2))^2-(nnM(j+1,3))^2);
            end
        end
    else
        for j=1:samples
            r=j/samples;
            c(1)=1;c(2)=-2*(a2+a3);c(3)=a2^2+4*a3*a2+a3^2-(b2^2+b3^2)/4/r^2;
            c(4)=-2*a2*a3*(a2+a3)+(a2*b3^2+a3*b2^2)/2/r^2;
            c(5)=a2^2*a3^2-(a2^2*b3^2+a3^2*b2^2)/4/r^2;
            rts=roots(c); realno=0;
            for k=1:4
                if abs(imag(rts(k)))<0.00001
                    realno=realno+1;
                end
            end
            rrts=zeros(realno,1);place=1;
            for k=1:4
                if abs(imag(rts(k)))<0.0001
                    rrts(place)=real(rts(k));
                    place=place+1;
                end
            end
            test=zeros(realno,1);
            for k=1:realno
                n=zeros(3,1);
                n(2)=b2/2/(rrts(k)-a2)/r;n(3)=b3/2/(rrts(k)-a3)/r;
                n=n/norm(n);
                test(k)=b2*n(2)+b3*n(3)+r*(-a1+a2*(n(2)^2-1)+a3*(n(3)^2-1));
            end
            [rdot(j+1,1),ind1]=max(test);
            [rdot(j+1,2),ind2]=min(test);
            nnM(j+1,2)=b2/2/(rrts(ind1)-a2)/r;nnM(j+1,3)=b3/2/(rrts(ind1)-a3)/r;
            nnm(j+1,2)=b2/(rrts(ind2)-a2)/r/2;nnm(j+1,3)=b3/2/(rrts(ind2)-a3)/r;
            if a1>a2
	    	fstar = b2^2/4/r/(a1-a2) + b3^2/4/r/(a1-a3) - r*(a2+a3);
            	if fstar > rdot(j+1,1) && r^2 >= b2^2/4/(a1-a2)^2 + b3^2/4/(a1-a3)^2
                	rdot(j+1,1) = fstar; 
                	nnM(j+1,2) = b2/2/r/(a1-a2);
                	nnM(j+1,3) = b3/2/r/(a1-a3);
                	nnM(j+1,1) = sqrt(1-(nnM(j+1,2))^2-(nnM(j+1,3))^2);
		end
            end
            nnM(j+1,1:3)=nnM(j+1,1:3)/norm(nnM(j+1,1:3));
	    nnm(j+1,1:3)=nnm(j+1,1:3)/norm(nnm(j+1,1:3));
            if rdot(j+1,1) <= 1E-14 && rdot(j,1) >= 0
               rcrit = (j - 1 + rdot(j,1)/(rdot(j,1)-rdot(j+1,1)))/samples; 
            end
        end
    end
    
% If all three b components are nonzero, we have to solve a two or four
% or six degree polynomial, depending on the multiplicity of the aj's 
    
elseif a1 == a3
    nnM(:,1)=b1*ones(samples+1,1)/bnorm;    nnm(:,1)=-b1*ones(samples+1,1)/bnorm;
    nnM(:,2)=b2*ones(samples+1,1)/bnorm;    nnm(:,2)=-b2*ones(samples+1,1)/bnorm;
    nnM(:,3)=b3*ones(samples+1,1)/bnorm;    nnm(:,3)=-b3*ones(samples+1,1)/bnorm;
    rdot(:,1) =  bnorm - 2*a1*(0:samples)/samples;
    rdot(:,2) = -bnorm - 2*a1*(0:samples)/samples;
    rcrit = bnorm /2/a1;
     
elseif xor(a1==a2,a2==a3)
    if a1==a2
        B1=sqrt(b1^2+b2^2);B2=b3;A1=a1;A2=a3;
    else
        B1=b1;B2=sqrt(b2^2+b3^2);A1=a1;A2=a3;
    end
    for j=1:samples
        r=j/samples;
        c(1)=1;c(2)=-2*(A1+A2);c(3)=A1^2+4*A1*A2+A2^2-(B1^2+B2^2)/4/r^2;
        c(4)=-2*A1*A2*(A1+A2)+(A1*B2^2+A2*B1^2)/2/r^2;
        c(5)=A1^2*A2^2-(A1^2*B2^2+A2^2*B1^2)/4/r^2;
        rts=roots(c);
        realno=0;
        for k=1:4
            if abs(imag(rts(k)))<0.00001
                realno=realno+1;
            end
        end
        rrts=zeros(realno,1);place=1;
        for k=1:4
            if abs(imag(rts(k)))<0.0001
                rrts(place)=real(rts(k));
           		place=place+1;
            end
        end
        test=zeros(realno,1);
        for k=1:realno
            n=zeros(3,1);
            n(1)=b1/(rrts(k)-a1)/2/r;n(2)=b2/2/(rrts(k)-a2)/r;n(3)=b3/2/(rrts(k)-a3)/r;
            n=n/norm(n);
            test(k)=b1*n(1)+b2*n(2)+b3*n(3)+r*(a1*(n(1)^2-1)+a2*(n(2)^2-1)+a3*(n(3)^2-1));
        end
        [rdot(j+1,1),ind1]=max(test);[rdot(j+1,2),ind2]=min(test);
        nnM(j,1)=b1/2/(rrts(ind1)-a1)/r;
	nnM(j,2)=b2/2/(rrts(ind1)-a2)/r;
	nnM(j,3)=b3/2/(rrts(ind1)-a3)/r;
	nnm(j,1)=b1/(rrts(ind2)-a1)/r/2;
	nnm(j,2)=b2/(rrts(ind2)-a2)/r/2;
	nnm(j,3)=b3/2/(rrts(ind2)-a3)/r;
        nnM(j,1:3)=nnM(j,1:3)/norm(nnM(j,1:3));
	nnm(j,1:3)=nnm(j,1:3)/norm(nnm(j,1:3));
        if rdot(j,1)>0&&rdot(j+1,1)<=1E-14
            rcrit=((j-1) - rdot(j,1)/(rdot(j+1,1)-rdot(j,1)))/samples;
        end
    end
    if rdot(samples,1)>=1
        rcrit=1;
    end       
else 
    c=zeros(7,1);
    for j=1:samples
        r=j/samples;
        c(1)=1;c(2)=-2*(a1+a2+a3);c(3)=a1^2+a2^2+a3^2+4*(a1*a2+a1*a3+a2*a3)-bnorm^2/4/r^2;
        c(4)=(b1^2*(a2+a3)+b2^2*(a3+a1)+b3^2*(a1+a2))/2/r^2-8*a1*a2*a3-...
		2*(a1*a2*(a1+a2)+a1*a3*(a1+a3)+a2*a3*(a2+a3));
        c(5)=a1^2*a2^2+a1^2*a3^2+a2^2*a3^2+4*a1*a2*a3*(a1+a2+a3)-(b1^2*(a2^2+4*a2*a3+a3^2)+...
		b2^2*(a1^2+4*a1*a3+a3^2)+b3^2*(a1^2+4*a1*a2+a2^2))/4/r^2;
        c(6)=(b1^2*a2*a3*(a2+a3)+b2^2*a1*a3*(a1+a3)+b3^2*a1*a2*(a1+a2))/2/r^2-...
		2*a1*a2*a3*(a1*a2+a1*a3+a2*a3);
        c(7)=a1^2*a2^2*a3^2-(b1^2*a2^2*a3^2+b2^2*a3^2*a1^2+b3^2*a1^2*a2^2)/4/r^2;
        rts=roots(c);
        realno=0;
        for k=1:6
            if abs(imag(rts(k)))<0.00001
                realno=realno+1;
            end
        end
        rrts=zeros(realno,1);place=1;
        for k=1:6
            if abs(imag(rts(k)))<0.0001
           		rrts(place)=real(rts(k));
            	place=place+1;
            end
        end
        test=zeros(realno,1);
        for k=1:realno
            n=zeros(3,1);
            n(1)=b1/(rrts(k)-a1)/2/r ;n(2)=b2/2/(rrts(k)-a2)/r ;n(3)=b3/2/(rrts(k)-a3)/r;
            n=n/norm(n);
            test(k)=b1*n(1)+b2*n(2)+b3*n(3)+r*(a1*(n(1)^2-1)+a2*(n(2)^2-1)+a3*(n(3)^2-1));
        end
        [rdot(j+1,1),ind1]=max(test);[rdot(j+1,2),ind2]=min(test);
        nnM(j+1,1)=b1/2/(rrts(ind1)-a1)/r;
	nnM(j+1,2)=b2/2/(rrts(ind1)-a2)/r;
	nnM(j+1,3)=b3/2/(rrts(ind1)-a3)/r;
	nnm(j+1,1)=b1/(rrts(ind2)-a1)/r/2;
	nnm(j+1,2)=b2/(rrts(ind2)-a2)/r/2;
	nnm(j+1,3)=b3/(rrts(ind2)-a3)/r/2;
    	nnM(j+1,:)=nnM(j+1,:)/norm(nnM(j+1,:));
        nnm(j+1,:)=nnm(j+1,:)/norm(nnm(j+1,:));
        if rdot(j,1)>0&&rdot(j+1,1)<=1E-14
            rcrit=(j - 1  - rdot(j,1)/(rdot(j+1,1)-rdot(j,1)))/samples;
        end
    end
    if rdot(samples+1,1)>=1
        rcrit=1;
    end     
end

% Time to plot.

if shouldplot == 1
    plot((0:samples)/samples,rdot(:,1));
    hold
    plot((0:samples)/samples,rdot(:,2));
    plot(0:1,[0,0],'Color','k');
    hold
    xlabel("r"); ylabel("dr/dt");
end

end
