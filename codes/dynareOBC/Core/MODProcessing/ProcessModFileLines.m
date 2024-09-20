function [ FileLines, Indices, StochSimulCommand, dynareOBC ] = ProcessModFileLines( FileLines, dynareOBC )
    dynareOBC.MaxFuncIndices = [];
    dynareOBC.NumberOfMax = 0;
    StochSimulCommand = '';
    Indices = struct;
    Indices.ModelStart = 0;
    Indices.ModelEnd = 0;
    Indices.ShocksStart = 0;
    Indices.ShocksEnd = 0;
    Indices.InitValStart = 0;
    Indices.InitValEnd = 0;
    Indices.SteadyStateModelStart = 0;
    Indices.SteadyStateModelEnd = 0;
    SearchState = 0;
    TempCounter = 0;
    i = 1;

    while i <= length( FileLines )
        write_i = i;
        line = FileLines{ i };
        switch SearchState
            case 0 % not in any block
                if ~isempty( regexp( line, '^model(\(.*\))?;', 'once' ) )
                    if Indices.ModelStart > 0
                        error( 'dynareOBC:DuplicateBlock', 'Duplicate model blocks detected.' );
                    end
                    line = regexprep( line, ',(use_dll|block|bytecode|linear)(?!\w)', '' );
                    line = regexprep( line, '(?<!\w)(use_dll|block|bytecode|linear),', '' );
                    line = regexprep( line, '\((use_dll|block|bytecode|linear)\)', '' );
                    FileLines{ write_i } = line;
                    SearchState = 1;
                    Indices.ModelStart = i;
                elseif ~isempty( regexp( line, '^shocks(\(.*\))?;', 'once' ) )
                    if Indices.ShocksStart > 0
                        error( 'dynareOBC:DuplicateBlock', 'Duplicate shocks blocks detected.' );
                    end
                    SearchState = 2;
                    Indices.ShocksStart = i;
                elseif ~isempty( regexp( line, '^initval(\(.*\))?;', 'once' ) )
                    if Indices.InitValStart > 0
                        error( 'dynareOBC:DuplicateBlock', 'Duplicate initval blocks detected.' );
                    end
                    SearchState = 3;
                    Indices.InitValStart = i;
                elseif ~isempty( regexp( line, '^steady_state_model(\(.*\))?;', 'once' ) )
                    if Indices.SteadyStateModelStart > 0
                        error( 'dynareOBC:DuplicateBlock', 'Duplicate steady_state_model blocks detected.' );
                    end
                    SearchState = 4;
                    Indices.SteadyStateModelStart = i;
                elseif ~isempty( regexp( line, '^stoch_simul(\(.*\))?(\s*\w+,?\s*)*;', 'once' ) )
                    StochSimulCommand = line;
                    PostScriptFileText = strjoin( FileLines( ( i + 1 ) : end ), '\n' );
                    PostScriptFile = fopen( 'dynareOBCTempPostScript.m', 'w' );
                    fprintf( PostScriptFile, '%s', PostScriptFileText );
                    fclose( PostScriptFile );
                    rehash;
                    FileLines = FileLines( 1:(i-1) );
                elseif ~isempty( regexp( line, '^var\>', 'once' ) ) || ~isempty( regexp( line, '^varexo\>', 'once' ) ) || ~isempty( regexp( line, '^varexo_det\>', 'once' ) ) || ~isempty( regexp( line, '^parameters\>', 'once' ) ) || ~isempty( regexp( line, '^trend_var\>', 'once' ) ) || ~isempty( regexp( line, '^log_trend_var\>', 'once' ) ) || ~isempty( regexp( line, '^model_local_variable\>', 'once' ) )
                    line = regexprep( line, '\$[^$]*\$', ' ' );
                    line = regexprep( line, '\(([^\)]*(''[^'']*'')?)*\)', ' ' );
                    line = regexprep( line, '\s+', ' ' );
                    FileLines{ write_i } = line;
                end
            case 1 % in the model block
                if strcmp( line, 'end;' )
                    SearchState = 0;
                    Indices.ModelEnd = i;
                else
                    line = regexprep( line, '\s*\[([^\]]*(''[^'']*'')?)*\]\s*', '' );
                    FileLines{ write_i } = line;
                    [ TempIndexStart, TempIndexEnd ] = regexp( line, '(?<=(^\#dynareOBCMaxFunc))\d+', 'once' );
                    if isempty( TempIndexStart )
                        [ FileLines, TempCounter, dynareOBC.NumberOfMax, write_i ] = ProcessModelLines( line, FileLines, TempCounter, dynareOBC.NumberOfMax, write_i );
                    else
                        dynareOBC.MaxFuncIndices( str2double( line( TempIndexStart:TempIndexEnd ) ) ) = i;
                    end
                end
            case 2 % in the shocks block
                if strcmp( line, 'end;' )
                    SearchState = 0;
                    Indices.ShocksEnd = i;
                end
            case 3 % in the initval block
                if strcmp( line, 'end;' )
                    SearchState = 0;
                    Indices.InitValEnd = i;
                end
            case 4 % in the steady_state_model block
                if strcmp( line, 'end;' )
                    SearchState = 0;
                    Indices.SteadyStateModelEnd = i;
                end
            otherwise
        end
        if write_i == i
            i = i + 1;
        end
    end
    if dynareOBC.Bypass
        dynareOBC.NumberOfMax = 0;
        dynareOBC.MaxFuncIndices = [];
    end

    if Indices.ModelStart == 0
        error( 'dynareOBC:MissingBlock', 'Start of model block was not found.' );
    end
    if Indices.ModelEnd == 0
        error( 'dynareOBC:MissingBlock', 'End of model block was not found.' );
    end
    if Indices.InitValStart > 0 && Indices.InitValEnd == 0
        error( 'dynareOBC:MissingBlock', 'End of initval block was not found.' );
    end
    if Indices.SteadyStateModelStart > 0 && Indices.SteadyStateModelEnd == 0
        error( 'dynareOBC:MissingBlock', 'End of steady_state_model block was not found.' );
    end
end
