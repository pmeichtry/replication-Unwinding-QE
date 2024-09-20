function [ x, M, oo ] = GlobalModelSolution( M, options, oo, dynareOBC )

    fprintf( '\n' );
    disp( 'Beginning to solve the fixed point problem.' );
    fprintf( '\n' );
    
    PI = dynareOBC.ParameterIndices_StateVariableAndShockCombinations(:);
    
    x = M.params( PI );
    if dynareOBC.Resume
        ResumeData = load( 'dynareOBCGlobalResume.mat' );
        ResumeParamNamesPI = ResumeData.ParamNamesPI;
        NewParamNames = cellstr( M.param_names );
        for i = 1 : length( ResumeParamNamesPI )
            ParamName = ResumeParamNamesPI{ i };
            j = find( strcmp( ParamName, NewParamNames ), 1 );
            if ~isempty( j )
                M.params( j ) = ResumeData.x( i );
            end
        end
        x = M.params( PI );
    end
    ParamNamesPI = cellstr( M.param_names( PI, : ) );
    
    try
        dynare_solve_path = which( 'dynare_solve.m' );
        dynare_solve_txt = fileread( dynare_solve_path );
        dynare_solve_txt = regexprep( dynare_solve_txt, 'options\.[dD]isplay\s*=\s*''iter''\s*;', 'options\.Display   =   ''off''   ;' );
        dynare_solve_file = fopen( dynare_solve_path, 'w' );
        fprintf( dynare_solve_file, '%s', dynare_solve_txt );
        fclose( dynare_solve_file );
        rehash;
    catch
        warning( 'dynareOBC:PatchDynareSolve', 'Error patching dynare_solve to disable output. We recommend you do this manually instead.' );
    end
    
    OpenPool;
    [ x, ResidualNorm ] = dynareOBC.FSolveFunctor( @(xx) GlobalModelSolutionInternal( xx, M, options, oo, dynareOBC ), x, 'OutputFcn', @( x, ~, ~ ) SaveProgress( x, ParamNamesPI ) );
    % x = cmaes( @(xx) norm( GlobalModelSolutionInternal( xx, M, options, oo, dynareOBC ) ), x, 2 );
    SaveProgress( x, ParamNamesPI );
    
    if ResidualNorm > sqrt( sqrt( eps ) )
        warning( 'dynareOBC:GlobalNonConvergence', 'Failed to find an exact fixed point. Returning a least squares deviation solution, which may compromise accuracy.' );
    end
        
    try
        dynare_solve_path = which( 'dynare_solve.m' );
        dynare_solve_txt = fileread( dynare_solve_path );
        dynare_solve_txt = regexprep( dynare_solve_txt, 'options\.[dD]isplay   =   ''off''   ;', 'options\.Display = ''iter'';' );
        dynare_solve_file = fopen( dynare_solve_path, 'w' );
        fprintf( dynare_solve_file, '%s', dynare_solve_txt );
        fclose( dynare_solve_file );
        rehash;
    catch
        warning( 'dynareOBC:PatchDynareSolve', 'Error patching dynare_solve to reenable output. We recommend you do this manually instead.' );
    end
    
    if any( ~isfinite( x ) )
        error( 'dynareOBC:GlobalBadParameters', 'Non-finite parameters were returned from the solution procedure.' );
    end
    
    M.params( PI ) = x;
    
    Info = -1;
    try
        [ dr, Info, M, ~, oo ] = resol( 0, M, options, oo );
        oo.dr = dr;
        oo.steady_state = oo.dr.ys;
    catch
    end
    
    if Info ~= 0
        error( 'dynareOBC:GlobalBadFinalPoint', 'Failed to solve the model at the final point.' );
    end
    
    x = reshape( x, [ size( dynareOBC.StateVariableAndShockCombinations, 1 ), dynareOBC.NumberOfMax ] );
end

function Stop = SaveProgress( x, ParamNamesPI )
    Stop = false;
    save dynareOBCGlobalResume.mat x ParamNamesPI;
end
