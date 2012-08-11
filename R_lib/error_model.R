#
#
# This sets up the R environment for the C calls
#
#
suppressPackageStartupMessages( library( mgcv ) )
dir.create( "error_model" )
logfname <- "error_model/r_output.log";
logfp = file( logfname, "w" );

MIN.NUM.OBS = 10;
DEFAULT.ERROR.EST = 0.001;
MIN.ERROR.PRB = 1e-6;

logit = function(x)log(x) - log(1-x)
logistic = function(x)exp(x)/(1+exp(x))

#
# Find the positions that we want to include in the model.
#
find.observed.predictors = function( predictor, cnts ) {
    agg = aggregate( cnts, by=list(predictor), FUN=sum );
    agg$Group.1[ which( agg$x > 0 ) ];
};

plot.hists = function( obs, pred, n=100 )
{
        p1 = hist( obs, plot=FALSE, n=n );
        p2 = hist( pred, plot=FALSE, n=n );        
        x_min = min(p1$mids, p2$mids);
        x_max = max(p1$mids, p2$mids);
        y_min = min(p1$density, p2$density);
        y_max = max(p1$density, p2$density);
        plot( p1$mids, p1$density, 
              xlim=c(x_min, x_max), ylim=c(y_min, y_max), 
              type='b', col='blue');
        points( p2$mids, p2$density, type='b', col='red');
        legend( x="topright", 
                legend=c("observed", "predicted"), 
                fill=c("blue", "red") );
        return();
}

build.mo = function( mm_cnts, cnts, pos, qual ) {
    # calculate bounds of the predictors. We do this to choose how complex to 
    # make the error model.
    obs.pos = find.observed.predictors(pos, cnts);
    obs.quals = find.observed.predictors(qual, cnts);

    # make sure that we have observed enough positions. It shouldn't be possible
    # for this not to be true because of the index probe size, so we throw a 
    # fatal exception if this doesn't hold.
    if( length( obs.pos ) < 12 ) {
        stop( "FATAL     : Too few observed unique positions (", 
              length(obs.pos), ") to build error model." );
        return;
    }

    #
    # Determine the model type
    #    
    response = cbind( mm_cnts, cnts-mm_cnts );
    if( length( obs.quals ) < 4 ) {
        eqn = formula( response ~ s(pos) );
        predictors = data.frame( pos=pos );
    }
    else if( FALSE && length( obs.quals ) < 6 ) {
        eqn = formula( response ~ s(pos) + qual*pos );
        predictors = data.frame( pos=pos );
    } else {
        eqn = formula( response ~ s(qual) ); # + s(pos) + s(qual, pos) );
        predictors = data.frame( pos=pos, qual=qual );
    }
    
    mo = gam( eqn, data=predictors, family=binomial );
            
    print( summary( mo ) );    
    
    mo;
}

# find a dense block subset of observations
find.data.subset = function( mm_cnts, cnts, pos, qual ) {
    obs.pos = find.observed.predictors(pos, cnts);
    valid.pos = ( pos >= min(obs.pos) & pos <= max(obs.pos) );
    obs.quals = find.observed.predictors(qual, cnts);
    valid.quals = ( qual >= min(obs.quals) & qual <= max(obs.quals) );
    
    # build the return values. Keep the boolean array to make it 
    # easier to rebuild the full matrices
    rv = list()
    rv$good.pos = valid.pos&valid.quals;
    rv$mm_cnts = mm_cnts[ valid.pos&valid.quals ];
    rv$cnts = cnts[ valid.pos&valid.quals ];
    rv$pos = pos[ valid.pos&valid.quals ];
    rv$qual = qual[ valid.pos&valid.quals ];
    
    rv;
}

predict_freqs_for_record = function( mm_cnts, cnts, pos, qual, plot.str=NULL ) {
    # set up the logging
    # sink( logfp, append=TRUE, type="message" );
    sink( logfp, append=TRUE );
    print( "Modeling error rates." );
    
    # If we havn't observed enough sequences  to do a good
    # job estimating the error rates, then use a default 
    # estimate ( but print out a warning ).
    n.obs = sum( cnts );
    print( paste( "Number of observations:", n.obs ) );
    if( n.obs < MIN.NUM.OBS )
    {
        print( paste( "WARNING     :  Number of observations (", n.obs, 
                      ") is too low. Using default error estimates.") );
        sink( NULL );
        return( mm_cnts/(cnts + 0.01) );         
    }
    
        
    # find a subset of the data to build the model on. This essentially
    # excludes predictors with 0 observations, but includes interior points
    data = find.data.subset( mm_cnts, cnts, pos, qual );
    
    mo = build.mo( data$mm_cnts, data$cnts, data$pos, data$qual );
    est.prbs = logistic( predict(mo) );
    
    rv = mm_cnts/(cnts + 0.01);
    pred.freqs = logistic( predict( mo ) );
    rv[] = max( pred.freqs );
    rv[ data$good.pos ] = pred.freqs;

    if( !is.null(plot.str) )
    {
        png(paste("error_model/", plot.str, ".png", sep=""), width=1200, height=1200);
        par( mfrow=c(2,2) );
        plot.hists( mm_cnts[cnts>0]/cnts[cnts>0], pred.freqs );
        qqplot( mm_cnts[cnts>0]/cnts[cnts>0], rv[cnts>0]);
        res = tryCatch({ 
            vis.gam( mo, plot.type="contour" );
        }, error = function(e) {
            hist(pred.freqs);
        })
        dev.off();
    }
    
    sink( NULL );
    
    return( rv );
}

predict_freqs = function( mm_cnts, cnts, pos, qual, plot.str=NULL ) {
     predict_freqs_for_record( mm_cnts, cnts, pos, qual, plot.str );
};