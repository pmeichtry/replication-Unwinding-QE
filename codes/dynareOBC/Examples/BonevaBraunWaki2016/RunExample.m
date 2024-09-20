addpath( '../FigureUtils' );

disp( 'This script illustrates IRFs to 5 standard deviation shocks in the Boneva Braun Waki (2016) model.' );

dynareOBC bbw2016.mod ShockScale=5

SaveFigure( [ ], 'ProductivityShock' );

h = groot;
h.CurrentFigure = h.Children( 2 );
SaveFigure( [ ], 'DemandShock' );

disp( 'Note that both shocks hit the zero lower bound.' );
disp( 'Note also from the output above that M is a P-matrix for this model, so there is a unique solution to the linearised model with the bound.' );
