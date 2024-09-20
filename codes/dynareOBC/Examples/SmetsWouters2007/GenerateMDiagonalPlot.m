figure( 100 );

subplot( 2, 2, 2 );

N = min( size( dynareOBC_.MMatrix ) );

plot( 0 : ( N - 1 ), diag( dynareOBC_.MMatrix ), 'k', 0 : ( N - 1 ), zeros( 1, N ), 'r' );

subplot( 2, 2, 1 );
YLim1 = get( gca, 'YLim' );
subplot( 2, 2, 2 );
YLim2 = get( gca, 'YLim' );

YLim = [ min( YLim1( 1 ), YLim2( 1 ) ), max( YLim1( 2 ), YLim2( 2 ) ) ];

subplot( 2, 2, 1 );
set( gca, 'YLim', YLim );
subplot( 2, 2, 2 );
set( gca, 'YLim', YLim );
