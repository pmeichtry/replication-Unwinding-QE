disp( 'We start by showing that M is not an S-matris in the Smets Wouters (2003) model, even with T=1000.' );

dynareOBC SW03.mod TimeToEscapeBounds=1000

disp( 'Observe that M was found to not be an S-matrix with T=1000.' );

disp( 'Things are different if we switch to a rule that responds to the price level and output growth.' );
disp( 'Press a key to continue:' );
pause;

dynareOBC SW03PLT.mod TimeToEscapeBounds=1000 FeasibilityTestGridSize=10 SkipQuickPCheck

disp( 'Observe that M was found to be an S-matrix for T=1000 and T=infinity. Observe also that M is a P-matrix for T=1000.' );
