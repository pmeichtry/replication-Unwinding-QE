function [ oo, dynareOBC ] = RunStochasticSimulation( M, options, oo, dynareOBC )

    % derived from simult.m
    PositiveVarianceShocks = setdiff( 1:dynareOBC.OriginalNumVarExo, find( diag(M.Sigma_e) < eps ) );
    NumberOfPositiveVarianceShocks = length( PositiveVarianceShocks );
    
    SqrtmSigma_e = spsqrtm( M.Sigma_e( PositiveVarianceShocks, PositiveVarianceShocks ) );

    if dynareOBC.SimulateOnGridPoints
        [U,D] = schur( full( dynareOBC.Var_z1 ), 'complex' );
        % assert( isreal( U ) );
        diagD = diag( D );
        % assert( isreal( diagD ) );
        RootD = sqrt( diagD );
        IDv = RootD > sqrt( eps );
        RootVar_z1 = U( :, IDv ) * diag( RootD( IDv ) );
        NumberOfPositiveVarianceVariables = size( RootVar_z1, 2 );
        QMCDraws = SobolSequence( NumberOfPositiveVarianceVariables + NumberOfPositiveVarianceShocks, dynareOBC.SimulationPeriods );
        ShockSequence = zeros( max( dynareOBC.OriginalNumVarExo, M.endo_nbr ), 2 * dynareOBC.SimulationPeriods );
        Z1Offsets = RootVar_z1 * QMCDraws( 1 : NumberOfPositiveVarianceVariables, : );
        ShockSequence( 1:M.endo_nbr, 1:2:end ) = Z1Offsets( oo.dr.inv_order_var, : );
        ShockSequence( PositiveVarianceShocks, 2:2:end ) = SqrtmSigma_e * QMCDraws( (NumberOfPositiveVarianceVariables+1):end, : );
        
        if ~isempty( dynareOBC.InitialStateFile )
            error( 'dynareOBC:SimulateOnGridPointsIncompatibleOptions', 'You cannot specify an initial state with SimulateOnGridPoints.' );
        end
        if ~isempty( dynareOBC.ShockSequenceFile )
            error( 'dynareOBC:SimulateOnGridPointsIncompatibleOptions', 'You cannot specify a sequence of shocks with SimulateOnGridPoints.' );
        end
    else
        if isempty( dynareOBC.ShockSequenceFile )
            ShockSequence = zeros( dynareOBC.OriginalNumVarExo, dynareOBC.SimulationPeriods );
            ShockSequence( PositiveVarianceShocks, : ) = SqrtmSigma_e * randn( NumberOfPositiveVarianceShocks, dynareOBC.SimulationPeriods );
        else
            if exist( dynareOBC.ShockSequenceFile, 'file' ) == 2
                FileData = load( dynareOBC.ShockSequenceFile, 'ShockSequence' );
                if isfield( FileData, 'ShockSequence' )
                    ShockSequence = FileData.ShockSequence;
                    if ~all( size( ShockSequence ) == [ dynareOBC.OriginalNumVarExo, dynareOBC.SimulationPeriods ] )
                        error( 'dynareOBC:LoadedShockSequenceWrongSize', 'The loaded ShockSequence was not the correct size. Expected %d x %d, found %d x %d.', dynareOBC.OriginalNumVarExo, dynareOBC.SimulationPeriods, size( ShockSequence, 1 ), size( ShockSequence, 2 ) );
                    end
                else
                    error( 'dynareOBC:LoadedShockSequenceFileWrongVariables', 'The given shock sequence file did not contain a ShockSequence variable.' );
                end
            else
                error( 'dynareOBC:NoShockSequenceFile', 'Failed to find the file: %s', dynareOBC.ShockSequenceFile );
            end
        end
    end
    
    if isempty( dynareOBC.InitialStateFile )
        Simulation = SimulateModel( ShockSequence, true, [], false, false );
    else
        nEndo = M.endo_nbr;
        
        if exist( dynareOBC.InitialStateFile, 'file' ) == 2
            FileData = load( dynareOBC.InitialStateFile, 'InitialState' );
            if isfield( FileData, 'InitialState' )
                InitialStateTmp = FileData.InitialState;
                if ~all( size( InitialStateTmp ) == [ nEndo - dynareOBC.NumberOfMax, 1 ] )
                    error( 'dynareOBC:LoadedInitialStateWrongSize', 'The loaded InitialState was not the correct size. Expected %d x %d, found %d x %d.', nEndo - dynareOBC.NumberOfMax, 1, size( InitialStateTmp, 1 ), size( InitialStateTmp, 2 ) );
                end
                InitialState = full( dynareOBC.Constant );
                InitialState( 1 : numel( InitialStateTmp ) ) = InitialStateTmp;
            else
                error( 'dynareOBC:LoadedInitialStateFileWrongVariables', 'The given shock sequence file did not contain an InitialState variable.' );
            end
        else
            error( 'dynareOBC:NoInitialStateFile', 'Failed to find the file: %s', dynareOBC.InitialStateFile );
        end

        nEndoZeros = zeros( nEndo, 1 );
        InitialFullState = struct;
        InitialFullState.bound_offset = nEndoZeros;
        InitialFullState.first = bsxfun( @minus, InitialState, full( dynareOBC.Constant ) );
        InitialFullState.total = InitialState;
        InitialFullState.total_with_bounds = InitialFullState.total;

        if dynareOBC.Order > 1
            InitialFullState.second = nEndoZeros;
            if dynareOBC.Order > 2
                InitialFullState.first_sigma_2 = nEndoZeros;
                InitialFullState.third = nEndoZeros;
            end
        end
        InitialFullState = orderfields( InitialFullState );
        
        Simulation = SimulateModel( ShockSequence, true, InitialFullState, false, false );
    end
    
    oo.exo_simul = ShockSequence';
    oo.endo_simul = Simulation.total_with_bounds;
    dynareOBC.SimulationsWithoutBounds = Simulation.total;
    
    EndoNames = strtrim( cellstr( M.endo_names ) );
    
    for i = 1 : M.endo_nbr
        assignin( 'base', EndoNames{ i }, oo.endo_simul( i, : ).' );
    end
    
    if dynareOBC.MLVSimulationMode > 0
        dynareOBC.MLVSimulationWithBounds = Simulation.MLVsWithBounds;
        dynareOBC.MLVSimulationWithoutBounds = Simulation.MLVsWithoutBounds;

        MLVNames = dynareOBC.MLVNames;
        MLVSelect = dynareOBC.MLVSelect;
        nMLVIRFs = length( MLVSelect );
        for i = 1 : nMLVIRFs
            MLVName = MLVNames{MLVSelect(i)};
            assignin( 'base', MLVName, Simulation.MLVsWithBounds.( MLVName )(:) );
        end
    end
    
    DispMoments( M, options, oo, dynareOBC );
    
end
