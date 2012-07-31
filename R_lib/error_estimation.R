get_beta_bounds = function( match, cnt ) {
    if( cnt == 0 ){
        return( c(0,1.0) );
    }
    
    lb = qbeta( 0.05, match+0.01, cnt-match );
    ub = qbeta( 0.95, match+0.01, cnt-match );
    
    c( lb, ub );
}

plot.marginal = function( t, index ) {    
    pos.mismatch.cnts = aggregate( t[[2]], by=list(t[[1]][,index]), FUN=sum );
    poss = pos.mismatch.cnts[[1]];
    pos.mismatch.cnts = pos.mismatch.cnts[[2]];
    pos.cnts = aggregate( t[[3]], by=list(t[[1]][,index]), FUN=sum )[[2]];
    bounds = mapply( get_beta_bounds, pos.mismatch.cnts, pos.cnts );

    plot( poss, pos.mismatch.cnts/pos.cnts, type='l', lwd=2, ylim=c( min(bounds), max(bounds) ) );
    lines( poss, t(bounds)[,1], type='l', lwd=0.5 );
    lines( poss, t(bounds)[,2], type='l', lwd=0.5 );        
}

plot.marginals = function( t ) {    
    par( mfrow=c(2,1) );
    plot.marginal( t, 1 );
    plot.marginal( t, 2 );
}
plot.marginals( t )

build_matrices_from_record = function( record ) {
    max_readlen = as.integer(record[4])
    min_qual    = as.integer(record[5])
    max_qual    = as.integer(record[6])
    nobs        = (max_qual-min_qual)*max_readlen
    
    mismatches = matrix( as.numeric(record[7:(nobs+6)]), ncol=max_readlen );
    cnts = matrix( as.numeric(record[(nobs+7):length(record)]), ncol=max_readlen );
    
    list( mismatches, cnts );
}

build_regression_vectors_from_record = function( record ) {
    max_readlen = as.integer(record[4]);
    min_qual    = as.integer(record[5]);
    max_qual    = as.integer(record[6]);
    nobs        = (max_qual-min_qual)*max_readlen;
    
    mismatches = as.numeric(record[7:(nobs+6)]);
    cnts = as.numeric(record[(nobs+7):length(record)]);    
    xs = data.frame( qual=rep(min_qual:(max_qual-1), max_readlen), 
                     pos=rep(1:max_readlen, times=max_qual-min_qual),
                     row.names=NULL );
    
    list( xs, mismatches, cnts );
}

build_nonzero_regression_vectors_from_record = function( record ) {
    data = build_regression_vectors_from_record( record );
    # find which types have zero observations
    zero_obs_indices = which( data[[3]] != 0 );
    list( data[[1]][zero_obs_indices,], 
          data[[2]][zero_obs_indices], 
          data[[3]][zero_obs_indices] );
}

predict = function( i, j, pos_mo, qual_mo, lambda ) {
    lambda*pos_mo[i] + (1-lambda)*qual_mo[j];
}

predict_all = function( pos_mo, qual_mo, lambda ) {
    rv = matrix(ncol=length(pos_mo), nrow=length(qual_mo))
    for( i in 1:length(pos_mo) ) 
    {
        for( j in 1:length(qual_mo) )
        {
            rv[j,i] = predict(i, j, pos_mo, qual_mo, lambda)
        }
    }
    rv;
}
predict_all( pos_mo$y, qual_mo$y, lambda )

est_loss = function( errors_mat, cnts_mat, pos_mo, qual_mo, lambda ) {
    total_loss = 0;
    freqs_mat = errors_mat/cnts_mat;
    for( i in 1:length(pos_mo) ) 
    {
        for( j in 1:length(qual_mo) )
        {
            if( !is.nan( freqs_mat[j,i] ) )
            {
                loss = sqrt((freqs_mat[j,i] - predict(i, j, pos_mo, qual_mo, lambda))**2);
                total_loss = total_loss + cnts_mat[j,i]*loss;
            }
        }
    }
    total_loss;
}

build_error_predictor = function( record ) {
    max_readlen = as.integer(record[4]);
    min_qualscore = as.integer(record[5]);
    max_qualscore = as.integer(record[6]);    

    mats = build_matrices_from_record( record)
    qual_freqs = colSums(t(mats[[1]]))/colSums(t(mats[[2]]))
    qual_cnts = colSums(t(mats[[2]]))
    pos_freqs = colSums(mats[[1]])/colSums(mats[[2]])
    pos_cnts = colSums(mats[[2]])
    
    pos_mo = smooth.spline( 1:length(pos_freqs), pos_freqs, w=pos_freqs*(1-pos_freqs) )
    qual_mo = smooth.spline( 1:length(qual_freqs), qual_freqs, w=qual_freqs*(1-qual_freqs) )
    
    lambdas = seq(0,1,0.01);
    lambda = lambdas[ which.min(sapply(lambdas, function(l) {
        est_loss(mats[[1]], mats[[2]], pos_mo$y, qual_mo$y, l )
    }))]
    
    par(mfrow=c(1,2));
    plot(qual_freqs);
    points( qual_mo, col='red' )
    plot(pos_freqs)
    points( pos_mo, col='red' )
    
    predict_error_rate = function( qualscore, position ) {
        position = min( max_readlen, position );
        qualscore = min( max_qualscore, qualscore );
        qualscore = max(1, qualscore-min_qualscore+1)
        predict( position, qualscore, pos_mo$y, qual_mo$y, lambda );
    }
    predict_error_rate;
}

if(0) {
    # try to do the spline computations by hand
    splineDesign(knots = 1:10, x = 4:7)
    
}

data = read.table( "error_stats.log", header=TRUE )
t = build_regression_vectors_from_record( data[1,] )
plot.marginals( t )



