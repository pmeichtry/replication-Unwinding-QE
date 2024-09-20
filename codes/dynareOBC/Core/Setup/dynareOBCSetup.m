function dynareOBCSetup( OriginalPath, CurrentFolder, dynareOBCPath, InputFileName, varargin )
%       This command runs dynareOBC with the model file specified in
%       the InputFileName arument.
%       Please type "dynareOBC help" to see the full instructions.
%
% INPUTS
%   InputFileName:  Input file name, "help", "addpath", "rmpath" or "testsolvers"
%   varargin:       List of arguments
%             
% OUTPUTS
%   none
%        
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2001-2015 Dynare Team and Tom Holden
%
% This file is part of dynareOBC.
%
% dynareOBC is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% dynareOBC is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with dynareOBC. If not, see <http://www.gnu.org/licenses/>.

    %% Check if need to download new files
    Update = true;
    
    CurrentDay = now;
    if exist( [ dynareOBCPath '/LastDependencyUpdate.mat' ], 'file' ) == 2
        LastUpdateStructure = load( [ dynareOBCPath '/LastDependencyUpdate.mat' ] );
        if isfield( LastUpdateStructure, 'CurrentDay' )
            if CurrentDay - LastUpdateStructure.CurrentDay < 1
                Update = false;
            end
        end
    end

    %% Initialization

    addpath( [ dynareOBCPath '/Core/' ] );
    addpath( [ dynareOBCPath '/Core/Utils/' ] );
    ClosePool( ismember( 'nopoolclose', lower( varargin ) ), true );
    dynareOBCCleanUp;
    
    addpath( [ dynareOBCPath '/Core/MChecks/' ] );
    addpath( [ dynareOBCPath '/Core/MODProcessing/' ] );
    
    addpath( [ dynareOBCPath '/Extern/nlma/' ] );
    addpath( [ dynareOBCPath '/Core/BaseSimulation/' ] );
    addpath( [ dynareOBCPath '/Core/Global/' ] );
    addpath( [ dynareOBCPath '/Core/ModelSolution/' ] );
    addpath( [ dynareOBCPath '/Core/OBCSimulation/' ] );
    addpath( [ dynareOBCPath '/Core/Display/' ] );
    
    addpath( [ dynareOBCPath '/Extern/EST-NLSS/' ] );
    addpath( [ dynareOBCPath '/Extern/DoubleDouble/' ] );
    addpath( [ dynareOBCPath '/Core/Estimation/' ] );
    
    addpath( [ dynareOBCPath '/Core/InnerProblem/' ] );
    
    DynareVersion = return_dynare_version( dynare_version );
    
    if DynareVersion < 4.4
        error( 'dynareOBC:OldDynare', 'Your version of dynare is too old to use with DynareOBC. Please update dynare.' );
    end
    
    if ~ismember( 'noclearall', lower( varargin ) )
        WarningState = warning( 'off', 'all' );
        try
            evalin( 'base', 'clear all;' );
        catch
        end
        try
            evalin( 'base', 'clear global;' );
        catch
        end
        try
            evalin( 'base', 'clearvars;' );
        catch
        end
        warning( WarningState );
    end
    
    DynarePath = fileparts( which( 'dynare' ) );
    NDynarePath = length( DynarePath );
    CurrentPaths = strsplit( path, ';' );
    for i = 1 : length( CurrentPaths )
        CurrentPath = CurrentPaths{i};
        if length( CurrentPath ) > NDynarePath && all( CurrentPath( 1:NDynarePath ) == DynarePath )
            rmpath( CurrentPath );
        end
    end
    addpath( DynarePath );

    CompileMEX( dynareOBCPath, Update );
    
    if strcmpi( InputFileName, 'addpath' )
        EnforceRequirementsAndGeneratePath( Update, OriginalPath, CurrentFolder, dynareOBCPath, InputFileName, varargin{:} );
        dynare_config;
        return;
    end

    global dynareOBC_;
    if isempty( dynareOBC_ )
        dynareOBC_ = struct;
    end
    
    dynareOBC_.DynareVersion = DynareVersion;

    FNameDots = strfind( InputFileName, '.' );
    if isempty( FNameDots )
        dynareOBC_.BaseFileName = InputFileName;
    else
        dynareOBC_.BaseFileName = InputFileName( 1:(FNameDots(end)-1) );
    end

    dynareOBC_ = SetDefaultOptions( dynareOBC_ );

    basevarargin = cell( 1, 0 );
    for i = 1:length( varargin )
        [ basevarargin, dynareOBC_ ] = ProcessArgument( varargin{ i }, basevarargin, dynareOBC_ );
    end

    if dynareOBC_.TimeToEscapeBounds <= 0
        error( 'dynareOBC:Arguments', 'TimeToEscapeBounds must be strictly positive.' );
    end
    if dynareOBC_.FirstOrderAroundRSS1OrMean2 > 2
        error( 'dynareOBC:Arguments', 'You cannot select both FirstOrderAroundRSS and FirstOrderAroundMean.' );
    end
    if dynareOBC_.QuasiMonteCarloLevel <= 0
        error( 'dynareOBC:Arguments', 'QuasiMonteCarloLevel must be strictly positive.' );
    end
    if dynareOBC_.CubatureRegions <= 0
        error( 'dynareOBC:Arguments', 'CubatureRegions must be strictly positive.' );
    end
    if dynareOBC_.CubatureCATCHDegree < 0
        error( 'dynareOBC:Arguments', 'CubatureCATCHDegree must be non-negative.' );
    end
    
    if dynareOBC_.FastCubature
        warning( 'dynareOBC:FastCubatureDeprecated', 'The FastCubature option has been deprecated. It is now equivalent to the Cubature option, which turns on cubature.' );
        dynareOBC_.Cubature = true;
    end
    
    if dynareOBC_.Global || dynareOBC_.CubatureRegions > 1 || dynareOBC_.CubatureCATCHDegree > 0 || dynareOBC_.HigherOrderSobolDegree > 0
        dynareOBC_.Cubature = true;
    end
    
    dynareOBC_ = rmfield( dynareOBC_, 'FastCubature' );
    dynareOBC_ = rmfield( dynareOBC_, 'ImportanceSampling' );
    dynareOBC_ = rmfield( dynareOBC_, 'NoCubature' );

    basevarargin( end + 1 : end + 5 ) = { 'noclearall', 'console', 'nograph', 'nointeractive', '-DdynareOBC=1' };
    
    if DynareVersion < 4.6
        basevarargin( end + 1 ) = { 'nolinemacro' };
    end

    if strcmpi( InputFileName, 'TestSolvers' )
        EnforceRequirementsAndGeneratePath( Update, OriginalPath, CurrentFolder, dynareOBCPath, InputFileName, varargin{:} );
        yalmiptest;
        if ~isempty( dynareOBC_.LPSolver )
            try
                yalmiptest( dynareOBC_.LPSolver );
            catch Error
                warning( 'dynareOBC:TestSolversError', '%s', Error.message );
            end
        end
        if ~isempty( dynareOBC_.MILPSolver )
            try
                yalmiptest( dynareOBC_.MILPSolver );
            catch Error
                warning( 'dynareOBC:TestSolversError', '%s', Error.message );
            end
        end
        Architecture = computer;
        if ( length( Architecture ) >= 5 ) && strcmp( Architecture(1:5), 'PCWIN' )
            try
                opti_Install_Test;
            catch Error
                warning( 'dynareOBC:TestSolversError', '%s', Error.message );
            end
        end
        return;
    end
    
    dynareOBC_.OtherMOD = [];
    
    CoreError = [];

    global M_ options_ oo_
    
    if ~isempty( dynareOBC_.OtherMODFile )
        
        if dynareOBC_.Cubature
            warning( 'dynareOBC:CubatureWithOtherMODFile', 'Enabling Cubature and OtherMODFile may produce unexpected behaviour. Agents will switch their long-run beliefs when there is a positive probability of no solution.' );
        end
        
        if dynareOBC_.Global
            dynareOBC_.Global = false;
            disp( ' ' );
            disp( 'Disabling Global since the OtherMODFile option is non-empty.' );
            disp( ' ' );
        end
        if dynareOBC_.Estimation
            dynareOBC_.Estimation = false;
            disp( ' ' );
            disp( 'Disabling Estimation since the OtherMODFile option is non-empty.' );
            disp( ' ' );
        end
        if dynareOBC_.Smoothing
            dynareOBC_.Smoothing = false;
            disp( ' ' );
            disp( 'Disabling Smoothing since the OtherMODFile option is non-empty.' );
            disp( ' ' );
        end
        if dynareOBC_.SimulateOnGridPoints
            dynareOBC_.SimulateOnGridPoints = false;
            disp( ' ' );
            disp( 'Disabling SimulateOnGridPoints since the OtherMODFile option is non-empty.' );
            disp( ' ' );
        end
        if dynareOBC_.SlowIRFs
            dynareOBC_.SlowIRFs = false;
            disp( ' ' );
            disp( 'Disabling SlowIRFs since the OtherMODFile option is non-empty.' );
            disp( ' ' );
        end
        if dynareOBC_.MedianIRFs
            dynareOBC_.MedianIRFs = false;
            disp( ' ' );
            disp( 'Disabling MedianIRFs since the OtherMODFile option is non-empty.' );
            disp( ' ' );
        end
        
        dynareOBCOriginal = dynareOBC_;
        
        dynareOBC_.SkipAllSimulation = true;
        
        OtherFileName = dynareOBC_.OtherMODFile;
        
        FNameDots = strfind( OtherFileName, '.' );
        if isempty( FNameDots )
            dynareOBC_.BaseFileName = OtherFileName;
        else
            dynareOBC_.BaseFileName = OtherFileName( 1:(FNameDots(end)-1) );
        end
        
        dynareOBC_.SaveMacroName = [ dynareOBC_.BaseFileName '-macroexp.mod' ];
        dynareOBC_.DataFile = '';
        
        dynareOBC_.OtherMODFile = InputFileName;
        
        OriginalOtherMODFileSwitchToProbability      = dynareOBC_.OtherMODFileSwitchToProbability;
        dynareOBC_.OtherMODFileSwitchToProbability   = dynareOBC_.OtherMODFileSwitchFromProbability;
        dynareOBC_.OtherMODFileSwitchFromProbability = OriginalOtherMODFileSwitchToProbability;

        try
            dynareOBC_ = dynareOBCCore( OtherFileName, basevarargin, dynareOBC_, @() EnforceRequirementsAndGeneratePath( Update, OriginalPath, CurrentFolder, dynareOBCPath, OtherFileName, varargin{:} ) );
        catch CaughtCoreError
            CoreError = CaughtCoreError;
        end
        
        dynareOBC_.SkipAllSimulation = false;
        
        Files = dir( '**/dynareOBCTemp*' );
        [ ~, FilesSortOrder ] = sort( cellfun( @length, { Files.folder } ), 'descend' );
        Files = Files( FilesSortOrder );
        for i = 1 : length( Files )
            File = Files( i );
            movefile( [ File.folder '/' File.name ], [ File.folder '/' strrep( File.name, 'dynareOBCTemp', 'dynareOBCOtherTemp' ) ], 'f' );
        end
        rehash path;
        
        dynareOBCOriginal.OtherMOD = struct;
        
        dynareOBCOriginal.OtherMOD.dynareOBC = dynareOBC_;
        dynareOBCOriginal.OtherMOD.M         = M_;
        dynareOBCOriginal.OtherMOD.options   = options_;
        dynareOBCOriginal.OtherMOD.oo        = oo_;
        
        dynareOBC_ = dynareOBCOriginal;
        
    end
    
    if isempty( CoreError )
    
        try
            dynareOBC_ = dynareOBCCore( InputFileName, basevarargin, dynareOBC_, @() EnforceRequirementsAndGeneratePath( Update, OriginalPath, CurrentFolder, dynareOBCPath, InputFileName, varargin{:} ) );
        catch CaughtCoreError
            CoreError = CaughtCoreError;
        end
    
    end
        
    
    %% Cleaning up

    if dynareOBC_.SaveMacro && ~isempty( dynareOBC_.SaveMacroName )
        try
            copyfile( 'dynareOBCTemp1.mod', dynareOBC_.SaveMacroName, 'f' );
        catch
        end
    end
    ClosePool( dynareOBC_.NoPoolClose, false );
    if ~dynareOBC_.NoCleanUp
        dynareOBCCleanUp;
    end

    evalin( 'base', 'global dynareOBC_' );
    
    if ~dynareOBC_.NoRestorePath
        path( OriginalPath );
    end
    
    if ~isempty( CoreError )
        rethrow( CoreError );
    end

end
