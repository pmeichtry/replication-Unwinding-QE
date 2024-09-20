function [ GlobalApproximationParameters, MaxArgValues, AmpValues ] = RunGlobalSolutionAlgorithm( basevarargin, SolveAlgo, FileLines, Indices, ToInsertBeforeModel, ToInsertInModelAtStart, ToInsertInModelAtEnd, ToInsertInInitVal, MaxArgValues, CurrentNumParams, CurrentNumVar, dynareOBC )
  
    MaxArgPattern = MaxArgValues( :, 1 ) > MaxArgValues( :, 2 );
    while true
        OldMaxArgValues = MaxArgValues;
        [ GlobalApproximationParameters, MaxArgValues, AmpValues ] = RunGlobalSolutionAlgorithmInternal( basevarargin, SolveAlgo, FileLines, Indices, ToInsertBeforeModel, ToInsertInModelAtStart, ToInsertInModelAtEnd, ToInsertInInitVal, MaxArgValues, MaxArgPattern, CurrentNumParams, CurrentNumVar, dynareOBC );
        
        NewMaxArgPattern = MaxArgValues( :, 1 ) > MaxArgValues( :, 2 );
        if any( NewMaxArgPattern ~= MaxArgPattern )
            MaxArgPattern = NewMaxArgPattern;
            MaxArgValues = OldMaxArgValues;
        else
            break;
        end
    end
    
end

function [ GlobalApproximationParameters, MaxArgValues, AmpValues ] = RunGlobalSolutionAlgorithmInternal( basevarargin, SolveAlgo, FileLines, Indices, ToInsertBeforeModel, ToInsertInModelAtStart, ToInsertInModelAtEnd, ToInsertInInitVal, MaxArgValues, MaxArgPattern, CurrentNumParams, CurrentNumVar, dynareOBC )

    fprintf( '\n' );
    disp( 'Generating the intermediate mod file.' );
    fprintf( '\n' );

    [ FileLines, ToInsertBeforeModel, ToInsertInModelAtEnd, ToInsertInInitVal, dynareOBC, AmpValues ] = ...
        InsertGlobalEquations( FileLines, ToInsertBeforeModel, ToInsertInModelAtEnd, ToInsertInInitVal, MaxArgValues, MaxArgPattern, CurrentNumParams, CurrentNumVar, dynareOBC );

    [ FileLines, Indices ] = PerformInsertion( ToInsertBeforeModel, Indices.ModelStart, FileLines, Indices );
    [ FileLines, Indices ] = PerformInsertion( ToInsertInModelAtStart, Indices.ModelStart + 1, FileLines, Indices );
    [ FileLines, Indices ] = PerformInsertion( ToInsertInModelAtEnd, Indices.ModelEnd, FileLines, Indices );
    [ FileLines, ~ ] = PerformInsertion( [ { 'initval;' } ToInsertInInitVal { 'end;' } ], Indices.ModelEnd + 1, FileLines, Indices );

    %Save the result

    FileText = strjoin( [ FileLines { [ 'stoch_simul(order=' int2str( dynareOBC.Order ) ',solve_algo=' int2str( SolveAlgo ) ',pruning,sylvester=fixed_point,irf=0,periods=0,nocorr,nofunctions,nomoments,nograph,nodisplay,noprint);' ] } ], '\n' ); % dr=cyclic_reduction,
    newmodfile = fopen( 'dynareOBCTempG.mod', 'w' );
    fprintf( newmodfile, '%s', FileText );
    fclose( newmodfile );

    fprintf( '\n' );
    disp( 'Calling dynare on the intermediate mod file.' );
    fprintf( '\n' );

    global M_ oo_ options_
    options_.solve_tolf = eps;
    options_.solve_tolx = eps;
    dynare( 'dynareOBCTempG.mod', basevarargin{:} );

    options_.solve_tolf = eps;
    options_.solve_tolx = eps;

    [ GlobalApproximationParameters, M, oo ] = GlobalModelSolution( M_, options_, oo_, dynareOBC );

    Generate_dynareOBCTempGetMaxArgValues( dynareOBC.DynareVersion, dynareOBC.NumberOfMax, 'dynareOBCTempG' );
    MaxArgValues = dynareOBCTempGetMaxArgValues( oo.steady_state, [ oo.exo_steady_state; oo.exo_det_steady_state ], M.params );
    if any( MaxArgValues( :, 1 ) == MaxArgValues( :, 2 ) )
        error( 'dynareOBC:JustBinding', 'DynareOBC does not support cases in which the constraint just binds in steady-state.' );
    end
        
end
