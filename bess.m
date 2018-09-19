function [] = bess(N)
% James Stine
% Oklahoma State University
% james.stine@okstate.edu
% Bessel Filter Root Extraction (Normalized)
% 
%
% The following meta file when executed in MATLAB, will output
% the coefficients for the B(s) for a Bessel Filter.  The
% polynomial is defined by
%
%               B(n) = [(2*n-1)*B(n-1)) + (s^2*B(n-2))]
%
%               B(0) = 1
%
%               B(1) = s + 1
%
% where n is the order of the Bessel Filter.  The algorithm basically
% performs two degree matrix shifting for the s^2 term, and adds the
% (2*n-1)*B(n-2) term.  Each term of the output represents the
% polynomial for each incremental order of s
%
% e.g.     (1 3 3) would be s^2 + 3*s + 3


a=[1];
b=[1 1];
% Changing u, changes the order of the Bessell Filter
for u=2:N

% This loop shifts the n-2 term of b two to the right
        for n=1:prod(size(a))
                for z = 1:2
                        h(1,z)=0;
                end
        h(1,n+z)=a(1,n);
        end


% This step multiplies the n-1 term of b by (2*n-1), where n is the order  
% of the Bessel Filter.  The order is defined as u above.
        e=b*(2*u-1);
        

% This step adds the the two processes state above.        
        for n=1:(prod(size(e)))
                h(1,n)=e(1,n)+h(1,n);
        end


% This step renames the n-1 term to the n-2 term for the next pass.  
% That is, if there is a next pass
        for n=1:(prod(size(b)))
                a(1,n)=b(1,n);
        end


% This step renames the current term to the n-1 term for the next pass.
% That is, if there is a next pass.
        for n=1:(prod(size(h)))
                b(1,n)=h(1,n);
        end
        
end

% This step flips the matrix 90 degress, because the algorithm above
% works the polynomials from right to left (that is, highest order
% from right to left).  If MATLAB is to perform the ROOTS function
% to find the roots of the polynomial, the polynomial must be from
% left to right (that is, highest order from right to left).
b=fliplr(h);

% This step finds the value to normalize our final result calculated
% above
K=(b(1,u))^(1/u)

% This step normalizes the roots
r=roots(b)/K;


rr=real(r)
ri=imag(r)


