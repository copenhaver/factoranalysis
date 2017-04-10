%%%
% A MATLAB implementation of Algorithm 1 from BCM17
% Written by Martin S. Copenhaver (www.mit.edu/~mcopen)
%%%

function [Theta,Phi] = FA(S,r,tol,maxiter)

    p = size(S,1);

    W = diag(rand(p,1));
    W = (p-r)*W/trace(W);
    
    iter = 1;
    pobj = Inf;
    cobj = 0;
    gap = Inf;
    
    while (iter < maxiter && gap > tol)

        cvx_begin sdp quiet
            variable Phi(p)
            minimize( trace(W*S)- trace(W*diag(Phi))  )
            subject to
                Phi >= 0;
                S - diag(Phi) >= 0;
        cvx_end
        
        A = S - diag(Phi);
        [U, SS] =svd(A);
        d=diag(U'*A*U);
        SS=SS*spdiags(sign(d.*diag(SS)),0,size(U,2),size(U,2));
        if ~(sort(diag(SS),'descend') == diag(SS))
            disp('Error - eigs not sorted!');
        end
        
        W = U*diag([zeros(r,1);ones(p-r,1)])*(U');
        
        iter = iter + 1;
        
        % update gap information
        
        pobj = cobj;
        cobj = trace(W*(S-diag(Phi)));
        gap = abs(cobj-pobj)/(abs(pobj)+.01)*100; % relative gap        
        
    end
    
    % compute Theta, which is closest rank-r matrix to S - diag(Phi)
    
    spec = diag(SS);
    Theta = U*diag([spec(1:r);zeros(p-r,1)])*U';

end
