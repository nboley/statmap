#####################################################
#
#  Normal spline regression 
#

test.spline_regression = function() {
    x = seq( 0, 10, 0.3 )
    y = sin( x ) + rnorm(0, 0.5, n=length(x) )

    basis = bs( x, knots=c(1,5,9) )
    # coefs = ((t(basis)%*%basis)**-1)%*%t(basis)%*%y
    # t(basis)%*%basis%*%Beta_hat = t(basis)%*%y
    coefs = solve( t(basis)%*%basis, t(basis)%*%y )
    y_hat = basis%*%coefs

    plot( y~x )
    points( y_hat~x, col='red', type='l' )
}


#####################################################
#
#  logistic spline regression 
#
library( splines );
library( gam )
logit = function(x)log(x) - log(1-x)
logistic = function(x)exp(x)/(1+exp(x))

fit.smooth.logistic.regression = function( x, cnts, n, B_old, knots ) {
    # build the basis function
    basis = bs( x, knots=knots );

    # estimate the current ps
    p_hat = logistic( basis%*%B_old );
    W = diag( as.vector(p_hat*(1-p_hat)) );
    
    ## Solve the update equation
    coefs = solve( t(basis)%*%W%*%basis, t(basis)%*%W%*%basis%*%B_old + t(basis)%*%(cnts/n-p_hat) );    
    coefs;
}

test.logistic.spline.regression = function() {
    x = seq( 0, 10, 0.3 )
    knots = c(2,5,7);
    
    # build a yector of probabilities
    prbs = sin( x ) + 1.01# + rnorm(0, 0.5, n=length(orig_x) )
    prbs[ which(prbs<0) ] = 0
    prbs = prbs/(max(prbs)*1.10)
    plot( x, prbs )
    
    # build a vactor of counts
    sample_size = 10
    cnts = sapply( prbs, function(p)rbinom(prob=p, size=sample_size, n=1) )
    plot( x, cnts )

    old_beta = rep( 0, 6 );
    for( loop in 1:100 )
    {
        beta = fit.smooth.logistic.regression( x, cnts, sample_size, old_beta, knots );
        if( sum((beta - old_beta)**2) < 1e-6 ) {
            print( paste( "NUM ITERATIONS: ", loop ) );
            break;
        }
        old_beta = beta;
    }

    # make a dense prediction
    new_x = seq(min(x),max(x),0.01);
    new_basis = bs( new_x, knots=knots );
    y_hat = sample_size*logistic(new_basis%*%beta);
    
    plot( x, cnts )
    points( x, sample_size*prbs, col='blue', type='l' )
    points( new_x, y_hat, type='l', col='red' )
    lines( x, 10*predict( gam( (y/10)~s(x, 6), link="logistic" ) ), col='green' )
}

test.logistic.spline.regression()





