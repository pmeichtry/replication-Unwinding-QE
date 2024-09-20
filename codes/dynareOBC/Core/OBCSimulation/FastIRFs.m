function [ oo, dynareOBC ] = FastIRFs( M, oo, dynareOBC )
% derived from nlma_irf.m
    
    IRFOffsets = struct;
    IRFsWithoutBounds = struct;
    T = dynareOBC.InternalIRFPeriods;
    Ts = dynareOBC.IRFPeriods;
    % Compute irf, allowing correlated shocks
    SS = M.Sigma_e + 1e-14 * eye( M.exo_nbr );
    cs = spsqrtm( SS );
    
    EndoNames = strtrim( cellstr( M.endo_names ) );
    ExoNames  = strtrim( cellstr( M.exo_names ) );
    
    if dynareOBC.MLVSimulationMode> 0
        VariableSelect = 1 : M.endo_nbr;
    else
        VariableSelect = dynareOBC.VariableSelect;
    end
    
    TempIRFLROffsets = repmat( dynareOBC.Mean, 1, T );
    TempIRFSROffsets = zeros( length( VariableSelect ), Ts );
    if dynareOBC.NumberOfMax > 0 && ( dynareOBC.Cubature || dynareOBC.Global )
        Shock = zeros( M.exo_nbr, 1 );
        y = FastIRFsInternal( Shock, 'the base path', M, oo, dynareOBC );
        for k = 1:length( VariableSelect )
            j = VariableSelect( k );
            TempIRFSROffsets( k, : ) = ( dynareOBC.MSubMatrices{ j }( 1:Ts, : ) * y )';
        end
    end
            
    for i = dynareOBC.ShockSelect
        Shock = zeros( M.exo_nbr, 1 ); % Pre-allocate and reset irf shock sequence
        Shock(:,1) = dynareOBC.ShockScale * cs( M.exo_names_orig_ord, i );
        
        [ y, TempIRFStruct ] = FastIRFsInternal( Shock, [ 'shock ' dynareOBC.Shocks{i} ], M, oo, dynareOBC );
        
        TempIRFs = TempIRFStruct.total - TempIRFLROffsets;
        
        for k = 1:length( VariableSelect )
            j = VariableSelect( k );
            IRFName = [ EndoNames{ j } '_' ExoNames{ i } ];
            CurrentIRF = TempIRFs( j, 1:Ts );
            IRFsWithoutBounds.( IRFName ) = CurrentIRF;
            if dynareOBC.NumberOfMax > 0
                CurrentIRF = CurrentIRF - TempIRFSROffsets( k, 1:Ts ) + ( dynareOBC.MSubMatrices{ j }( 1:Ts, : ) * y )';
            end
            oo.irfs.( IRFName ) = CurrentIRF;
            IRFOffsets.( IRFName ) = TempIRFLROffsets( j, 1:Ts );
            assignin( 'base', IRFName, CurrentIRF.' );
        end
    end
    dynareOBC.IRFOffsets = IRFOffsets;
    dynareOBC.IRFsWithoutBounds = IRFsWithoutBounds;

    if dynareOBC.MLVSimulationMode > 0
        
        if dynareOBC.Cubature
            fprintf( '\n' );
            disp( 'MLVSimulationMode>0 is not supported with cubature and without the SlowIRFs option. Skipping MLV simulation.' );
            disp( 'Consider specifying SlowIRFs in future.' );
            fprintf( '\n' );
            return
        end
        
        if dynareOBC.Global
            fprintf( '\n' );
            disp( 'MLVSimulationMode>0 is not supported with global and without the SlowIRFs option. Skipping MLV simulation.' );
            disp( 'Consider specifying SlowIRFs in future.' );
            fprintf( '\n' );
            return
        end
        
        if dynareOBC.MLVSimulationMode > 1
            fprintf( '\n' );
            disp( 'Only MLVSimulationMode=1 is supported without the SlowIRFs option. DynareOBC will act as if you specified MLVSimulationMode=1.' );
            disp( 'Consider specifying SlowIRFs in future.' );
            fprintf( '\n' );            
        end
        
        if dynareOBC.Order > 1
            fprintf( '\n' );
            disp( 'MLV simulation without SlowIRFs is a bad idea with order>1.' );
            disp( 'Consider specifying SlowIRFs in future.' );
            fprintf( '\n' );            
        end
        
        MLVNames = dynareOBC.MLVNames;
        nMLV = length( MLVNames );
        OriginalVarSelect = false( M.endo_nbr, 1 );
        OriginalVarSelect( 1:dynareOBC.OriginalNumVar ) = true;
        
        ParamVec = M.params;
        SteadyState = full( oo.dr.ys( 1:dynareOBC.OriginalNumVar ) );

        if dynareOBC.OriginalMaximumEndoLag > 0
            LagIndices = dynareOBC.OriginalLeadLagIncidence( dynareOBC.OriginalMaximumEndoLag, : ) > 0;
        else
            LagIndices = [];
        end
        CurrentIndices = dynareOBC.OriginalLeadLagIncidence( dynareOBC.OriginalMaximumEndoLag + 1, : ) > 0;
        if size( dynareOBC.OriginalLeadLagIncidence, 1 ) >= dynareOBC.OriginalMaximumEndoLag + 2
            LeadIndices = dynareOBC.OriginalLeadLagIncidence( dynareOBC.OriginalMaximumEndoLag + 2, : ) > 0;
        else
            LeadIndices = [];
        end
        
        FutureValues = nan( sum( LeadIndices ), 1 );

        LagValuesWithBounds = dynareOBC.Mean( OriginalVarSelect );
        LagValuesWithBoundsLagIndices = LagValuesWithBounds( LagIndices );
        CurrentValuesWithBoundsCurrentIndices = LagValuesWithBounds( CurrentIndices );

        MLVValuesWithBounds = dynareOBCTempGetMLVs( full( [ LagValuesWithBoundsLagIndices; CurrentValuesWithBoundsCurrentIndices; FutureValues ] ), zeros( size( Shock ) ), ParamVec, SteadyState );

        TempIRFLROffsets = zeros( nMLV, Ts );
        for j = 1 : nMLV
            TempIRFLROffsets( j, : ) = MLVValuesWithBounds( j );
        end

        for i = dynareOBC.ShockSelect
            
            ShockName = ExoNames{ i };
            
            IRFsAsArrayWithBounds = zeros( dynareOBC.OriginalNumVar, Ts );
            IRFsAsArrayWithoutBounds = zeros( dynareOBC.OriginalNumVar, Ts );
            
            for j = 1 : dynareOBC.OriginalNumVar
                IRFName = [ EndoNames{ j } '_' ShockName ];
                IRFsAsArrayWithBounds( j, : ) = oo.irfs.( IRFName ) + IRFOffsets.( IRFName );
                IRFsAsArrayWithoutBounds( j, : ) = IRFsWithoutBounds.( IRFName ) + IRFOffsets.( IRFName );
            end
        
            MLVsWithBounds = struct;
            MLVsWithoutBounds = struct;
            LagValuesWithBounds = dynareOBC.Mean( OriginalVarSelect );
            LagValuesWithoutBounds = LagValuesWithBounds;

            p = TimedProgressBar( ceil( Ts / 10 ), 50, [ 'Computing model local variable paths for response to shock ' ShockName '. Predicted to finish within ' ], '. Progress: ', [ 'Computing model local variable paths for response to shock ' ShockName '.               Completed in ' ] );

            for j = 1 : nMLV
                MLVName = MLVNames{j};
                MLVsWithBounds.( MLVName ) = NaN( 1, Ts );
                MLVsWithoutBounds.( MLVName ) = NaN( 1, Ts );
            end

            for t = 1 : Ts
                % clear the last warning
                lastwarn( '' );
                % temporarily disable warnings
                WarningState = warning( 'off', 'all' );
                % wrap in a try catch block to ensure they're re-enabled
                try
                    CurrentValuesWithBounds = IRFsAsArrayWithBounds( :, t );
                    CurrentValuesWithoutBounds = IRFsAsArrayWithoutBounds( :, t );

                    if t == 1
                        CurrentShock = Shock;
                    else
                        CurrentShock = zeros( size( Shock ) );
                    end
                    
                    LagValuesWithBoundsLagIndices = LagValuesWithBounds( LagIndices );
                    LagValuesWithoutBoundsLagIndices = LagValuesWithoutBounds( LagIndices );
                    CurrentValuesWithBoundsCurrentIndices = CurrentValuesWithBounds( CurrentIndices );
                    CurrentValuesWithoutBoundsCurrentIndices = CurrentValuesWithBounds( CurrentIndices );

                    MLVValuesWithBounds = dynareOBCTempGetMLVs( full( [ LagValuesWithBoundsLagIndices; CurrentValuesWithBoundsCurrentIndices; FutureValues ] ), CurrentShock, ParamVec, SteadyState );
                    if dynareOBC.NumberOfMax > 0
                        MLVValuesWithoutBounds = dynareOBCTempGetMLVs( full( [ LagValuesWithoutBoundsLagIndices; CurrentValuesWithoutBoundsCurrentIndices; FutureValues ] ), CurrentShock, ParamVec, SteadyState );
                    else
                        MLVValuesWithoutBounds = MLVValuesWithBounds;
                    end

                    for j = 1 : nMLV
                        MLVName = MLVNames{j};
                        MLVsWithBounds.( MLVName )( t ) = MLVValuesWithBounds( j );
                        MLVsWithoutBounds.( MLVName )( t ) = MLVValuesWithoutBounds( j );
                    end
                    
                    LagValuesWithBounds = CurrentValuesWithBounds;
                    LagValuesWithoutBounds = CurrentValuesWithoutBounds;
                catch Error
                    warning( WarningState );
                    rethrow( Error );
                end

                if ~isempty( p ) && rem( t, 10 ) == 0
                    p.progress;
                end
            end
            if ~isempty( p )
                p.stop;
            end
            
            for j = 1 : nMLV
                MLVName = MLVNames{j};
                IRFName = [ MLVName '_' ShockName ];
                oo.irfs.( IRFName ) = MLVsWithBounds.( MLVName ) - TempIRFLROffsets( j, 1:Ts );
                IRFsWithoutBounds.( IRFName ) = MLVsWithoutBounds.( MLVName ) - TempIRFLROffsets( j, 1:Ts );
                IRFOffsets.( IRFName ) = TempIRFLROffsets( j, 1:Ts );
                assignin( 'base', IRFName, oo.irfs.( IRFName ).' );
            end
        end
        dynareOBC.IRFOffsets = IRFOffsets;
        dynareOBC.IRFsWithoutBounds = IRFsWithoutBounds;
    end
    
end

function [ y, TempIRFStruct ] = FastIRFsInternal( Shock, ShockName, M, oo, dynareOBC )
    TempIRFStruct = ExpectedReturn( Shock, M, oo.dr, dynareOBC );

    UnconstrainedReturnPath = vec( TempIRFStruct.total( dynareOBC.VarIndices_ZeroLowerBounded, : )' );
    
    if dynareOBC.NumberOfMax > 0
        if dynareOBC.Cubature
            [ y, GlobalVarianceShare ] = PerformCubature( UnconstrainedReturnPath, oo, dynareOBC, TempIRFStruct.first, false, [ 'Computing required integral for fast IRFs for ' ShockName '. Predicted to finish within ' ], '. Progress: ', [ 'Computing required integral for fast IRFs for ' ShockName '.               Completed in ' ] );
            if dynareOBC.Global
                y = SolveGlobalBoundsProblem( y, GlobalVarianceShare, UnconstrainedReturnPath, TempIRFStruct.total( dynareOBC.VarIndices_ZeroLowerBoundedLongRun, : )', dynareOBC );
            end
        else
            y = SolveBoundsProblem( UnconstrainedReturnPath );
            if ~isempty( dynareOBC.IRFsForceAtBoundIndices ) || ~isempty( dynareOBC.IRFsForceNotAtBoundIndices )
                seps = sqrt( eps );
                y( dynareOBC.IRFsForceNotAtBoundIndices ) = 0;
                Periods = union( find( y ~= 0 ), dynareOBC.IRFsForceAtBoundIndices );
                y( Periods ) = -dynareOBC.MMatrix( Periods, Periods ) \ UnconstrainedReturnPath( Periods );
                while any( UnconstrainedReturnPath( dynareOBC.IRFsForceAtBoundIndices ) + dynareOBC.MMatrix( dynareOBC.IRFsForceAtBoundIndices, : ) * y > seps ) || any( UnconstrainedReturnPath + dynareOBC.MMatrix * y < -seps )
                    OldPeriods = Periods;
                    Periods = setdiff( union( Periods, find( UnconstrainedReturnPath( dynareOBC.sIndices ) + dynareOBC.MsMatrix * y < seps ) ), dynareOBC.IRFsForceNotAtBoundIndices );
                    if ( numel( Periods ) == numel( OldPeriods ) ) && all( Periods(:) == OldPeriods(:) )
                        error( 'dynareOBC:ForceAtBoundsIterativeFailure', 'The iterative ForceAtBounds/ForceNotAtBounds procedure failed. Note, this does not mean very much about the model!' );
                    end
                    y( Periods ) = -dynareOBC.MMatrix( Periods, Periods ) \ UnconstrainedReturnPath( Periods );
                end
                if dynareOBC.DisplayBoundsSolutionProgress
                    disp( y );
                end
            end
        end
    else
        y = [];
    end
end
