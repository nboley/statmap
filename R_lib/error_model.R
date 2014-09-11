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


predict.freqs = function( mm.cnts, cnts, pos, qual, do.plot=FALSE ) {
      aggregate.marginal.data = function( mm.cnts, cnts, by ) {
        mm.cnts.agg = aggregate( mm.cnts, by=list( by ), sum );
        cnts.agg = aggregate( cnts, by=list( by ), sum );
        data.frame(mm.cnts=mm.cnts.agg$x, cnts=cnts.agg$x, pred=cnts.agg$Group.1);
      };
      
      initial.smoothing = function( mm.cnts, cnts ) {
        # if we don't have enough observations, then
        # use logistic regression ( since we only have 1 or 2 points
        # the linear model is as good as we can do.
        if( sum( cnts > 0 ) == 1 )
        {
            return( rep( (mm.cnts[cnts>0]+1)/(cnts[cnts>0]+40), length(cnts) ) );
        } else if( sum(cnts > 0) < 3 )
        {
          response = cbind( mm.cnts+1, cnts-mm.cnts+20 );
          x = 1:length(cnts);
          mo = glm( response~x, family=binomial )
          print( summary( mo ) );
          return( predict( mo ) );
        }
        prbs = mm.cnts/(cnts + 1e-6);
        prbs = sapply( prbs, function(m)max( MIN.ERROR.PRB, m ) );

        # var(p) = np(1-p) => var(logit(p)) ~= n
        weights = sqrt( cnts );
        mo = smooth.spline( logit(prbs), w=weights );
        predict( mo )$y;
      };
      
      build.marginal = function(mm.cnts, cnts, pred) {
        agg = aggregate.marginal.data( mm.cnts, cnts, pred );
        pred.smooth = data.frame( pred=agg$pred,
          pred.prb=initial.smoothing( agg$mm, agg$cnts ));
        
        if( do.plot ) {
          obs.log.prbs = log10((agg$mm+0.0001)/(agg$cnts+0.01))
          plot( agg$pred, obs.log.prbs, xlab="Position in Read", ylab="log10 Predicted Mismatch Rate" );
          lines( agg$pred, log10(logistic(pred.smooth$pred.prb)), col='red' );
        }

        return( pred.smooth );
      }
      
      data = data.frame(mm.cnts=mm.cnts, cnts=cnts, pos=pos, qual=qual);
      pos.smooth.marginal = build.marginal(mm.cnts, cnts, pos)
      qual.smooth.marginal = build.marginal(mm.cnts, cnts, qual)
      
      data = merge( data, qual.smooth.marginal, by.x="qual", by.y="pred" )
      names(data)[5] <- "qual.smooth.prb";
      
      data = merge( data, pos.smooth.marginal, by.x="pos", by.y="pred" )
      names(data)[6] <- "pos.smooth.prb";
      
      response = cbind( data$mm.cnts, data$cnts-data$mm.cnts );
      if( nrow(qual.smooth.marginal) > 2 ) {
        error.formula = formula( response~qual.smooth.prb*pos.smooth.prb );
      } else {
        error.formula = formula( response~pos.smooth.prb );
      }
      mo = glm( error.formula, data=data, family=binomial );
      print( summary( mo ) );
      
      data$pred.freqs = logistic( predict( mo ) );
      data;
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
    
    new.data = predict.freqs( data$mm.cnts, data$cnts,
      data$pos, data$qual, !is.null(plot.str) );

    min.error.rate = max(MIN.ERROR.PRB, min(new.data$pred.freqs));
    rv = rep( min.error.rate, length(data$good.pos) );
    rv[ data$good.pos ] = new.data$pred.freqs;
    rv[ rv < min.error.rate ] = min.error.rate;
    
    if( !is.null(plot.str) )
    {
      # set the seed so that Identical models yield identical plots
      #set.seed( 0 );
      plot( ecdf( new.data$mm.cnts/new.data$cnts ), xlab="Reads Mismatch Fraction", col='blue' );
      for( i in 1:10 ) {
        r.sample = rbinom( length(new.data$cnts), new.data$cnts, new.data$pred.freqs );
        lines( ecdf( r.sample/new.data$cnts ), col='red' );      
      }
      lines( ecdf( new.data$mm.cnts/new.data$cnts ), col='blue' );
      
      qqplot( new.data$mm.cnts/new.data$cnts, r.sample/new.data$cnts, xlab="Observed Mismatch Fraction", ylab="Predicted Mismatch Fraction" );
      abline( a=0, b=1 )
        
      dev.off();
    }
    
    sink( NULL );
    
    return( rv );
}

predict_freqs = function( mm.cnts, cnts, pos, qual, plot.str=NULL ) {
     predict_freqs_for_record( mm.cnts, cnts, pos, qual, plot.str );
};
