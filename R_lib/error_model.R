#
#
# This sets up the R environment for the C calls
#
#
dir.create( "error_model" )
logfname <- "error_model/r_output.log";
logfp = file( logfname, "w" );

MIN.NUM.OBS = 10;
MIN.ERROR.PRB = 1e-4;

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

build.mo = function( mm.cnts, cnts, pos, qual, do.plot=FALSE ) {
      aggregate.marginal.data = function( mm.cnts, cnts, by ) {
        mm.cnts.agg = aggregate( mm.cnts, by=list( by ), sum );
        cnts.agg = aggregate( cnts, by=list( by ), sum );
        data.frame(mm.cnts=mm.cnts.agg$x, cnts=cnts.agg$x, pred=cnts.agg$Group.1);
      };
      
      initial.smoothing = function( mm.cnts, cnts ) {
        prbs = mm.cnts/(cnts + 1e-6);
        prbs = sapply( prbs, function(m)max( MIN.ERROR.PRB, m ) );
        # weights = sqrt( prbs*(1-prbs)*cnts );
        weights = sqrt(cnts*prbs);
        mo = smooth.spline( prbs, w=weights );
        predict( mo )$y;
      };
      
      build.marginals = function(data) {
        pos.agg = aggregate.marginal.data( data$mm.cnts, data$cnts, data$pos );
        pos.smooth = data.frame( pos=pos.agg$pred,
          pos.prb=initial.smoothing( pos.agg$mm, pos.agg$cnts ));

        if( do.plot ) {
          plot( pos.agg$pred, pos.agg$mm/(pos.agg$cnts+0.01) );
          lines( pos.agg$pred, pos.smooth$pos.prb, col='red' );
        }
        
        qual.agg = aggregate.marginal.data( data$mm.cnts, data$cnts, data$qual );
        qual.smooth = data.frame( qual=qual.agg$pred,
          qual.prb=initial.smoothing( qual.agg$mm, qual.agg$cnts ));
        
        if( do.plot ) {
          plot( qual.agg$pred, qual.agg$mm/(qual.agg$cnts+0.01) );
          lines( qual.agg$pred, qual.smooth$qual.prb, col='red' );
        }

        data.m1 = merge( data, pos.smooth );
        data.m2 = merge( data.m1, qual.smooth );
        return( data.m2 );
      }
      
      data = data.frame(mm.cnts=mm.cnts, cnts=cnts, pos=pos, qual=qual);
      new.data = build.marginals( data )
      print( head( new.data ) );
      
      response = cbind( new.data$mm.cnts, new.data$cnts-new.data$mm.cnts );
      mo = glm( response~qual.prb, data=new.data, family=binomial );
      
      print( summary( mo ) );

      mo;
}

# find a dense block subset of observations
find.data.subset = function( mm.cnts, cnts, pos, qual ) {
    obs.pos = find.observed.predictors(pos, cnts);
    valid.pos = ( pos >= min(obs.pos) & pos <= max(obs.pos) );
    obs.quals = find.observed.predictors(qual, cnts);
    valid.quals = ( qual >= min(obs.quals) & qual <= max(obs.quals) );
    
    # build the return values. Keep the boolean array to make it 
    # easier to rebuild the full matrices
    rv = list()
    rv$good.pos = valid.pos&valid.quals;
    rv$mm.cnts = mm.cnts[ valid.pos&valid.quals ];
    rv$cnts = cnts[ valid.pos&valid.quals ];
    rv$pos = pos[ valid.pos&valid.quals ];
    rv$qual = qual[ valid.pos&valid.quals ];
    
    rv;
}

predict_freqs_for_record = function( mm.cnts, cnts, pos, qual, plot.str=NULL ) {
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
        return( mm.cnts/(cnts + 0.01) );         
    }

    # find a subset of the data to build the model on. This essentially
    # excludes predictors with 0 observations, but includes interior points
    data = find.data.subset( mm.cnts, cnts, pos, qual );

    # Set up the plotting environment
    if( !is.null(plot.str) )
    {
      save( data, file=paste("error_model/", plot.str, ".obj", sep="") );

      plot.fname = paste("error_model/", plot.str, ".png", sep="");
      png(plot.fname, width=1200, height=1200);
      par( mfrow=c(2,2) );
    }    
    
    mo = build.mo( data$mm.cnts, data$cnts,
      data$pos, data$qual, !is.null(plot.str) );
    
    pred.freqs = logistic( predict( mo ) );
    rv = rep( max(MIN.ERROR.PRB, min(pred.freqs)), length(data$good.pos) );
    rv[ data$good.pos ] = pred.freqs;
    
    if( !is.null(plot.str) )
    {
      # set the seed so that Identical models yield identical plots
      set.seed( 0 );
      r.sample = rbinom( length(cnts), cnts, rv );
      plot.hists( mm.cnts[cnts>0]/cnts[cnts>0], r.sample[cnts>0]/cnts[cnts>0] );
      qqplot( mm.cnts[cnts>0]/cnts[cnts>0], r.sample[cnts>0]/cnts[cnts>0]);
      abline( a=0, b=1 )
        
      dev.off();
    }
    
    sink( NULL );
    
    return( rv );
}

predict_freqs = function( mm.cnts, cnts, pos, qual, plot.str=NULL ) {
     predict_freqs_for_record( mm.cnts, cnts, pos, qual, plot.str );
};
