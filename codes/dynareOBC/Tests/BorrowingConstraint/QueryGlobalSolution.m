function [ Xnew, DXnew, DDXnew ] = QueryGlobalSolution( B, A )
    persistent V X PP Bv Av beta mu rho sigma Ybar R
    
    if coder.target('MATLAB')
        if isempty( Bv )
            Results = load( 'GlobalResults.mat' );
            V = Results.V;
            X = Results.X;
            Bv = Results.Bv;
            PP = pchip( Bv, V ); % calculate rather than using saved as MATLAB Coder's version of pchip produces a different structure
            Av = Results.Av;
            beta = Results.beta;
            mu = Results.mu;
            rho = Results.rho;
            sigma = Results.sigma;
            Ybar = Results.Ybar;
            R = Results.R;
        end
    else
        Results = coder.load( 'GlobalResults.mat' );
        V = Results.V;
        X = Results.X;
        PP = Results.PP;
        Bv = Results.Bv;
        Av = Results.Av;
        beta = Results.beta;
        mu = Results.mu;
        rho = Results.rho;
        sigma = Results.sigma;
        Ybar = Results.Ybar;
        R = Results.R;
    end
    [ ~, Xnew ] = EvaluateValueFunctionOffGrid( B, A, Bv, Av, PP, max( max( V ) ), X, beta, Ybar, R, mu, rho, sigma );
    DXnew = NaN( 2, 1 );
    DDXnew = NaN( 2, 2 );
end
