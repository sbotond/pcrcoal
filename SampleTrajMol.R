##  
## Copyright 2010- Botond Sipos  
## See the package description for licensing information.   
## 

# Sample a size trajectories for each sub-population 
# (descendants of a single starting molecule)
.sampleTrajectory<-function(this){

# Check if everything is sane about the object:
.checkInput(this)

# Construct trajectory matrix: 
        traj.mat<-matrix(
            nrow=this@initial.size,	# a row for each starting molecule
            ncol=this@nr.cycles + 1  	# a column for each "generation"
        );

tries  <- 1

# Try to sample a sane subsample size vector. Give up after max.tries times.
while (TRUE) {

    if(tries >= this@max.tries){
        stop("\n\nTried ",tries," times to sample trajectories, but still could not obtain a large enough final pool size.\nAborting simulation. You should try again with higher efficiencies/cycle numbers/max.tries.\n\n")
    }

    # Initialize molecules:
    traj.mat[,1]<-1;

    for(i in 1:nrow(traj.mat)){
        for(j in 2:ncol(traj.mat)){
        
            prev.size <- traj.mat[i, j-1]
            traj.mat[i,j]<- prev.size + rbinom(n=1, size=prev.size,prob=this@efficiencies[j-1]) 
	    print(this@efficiencies[j-1])

        }

    }

    # Sanity checks:
    if(any(is.nan(traj.mat))){
        stop("NaN elements in the trajectory matrix! Maybe the cycle number is too high?")
    }

    # Check if final size is large enough:
    final.size<-sum(traj.mat[,ncol(traj.mat)])
    if(this@sample.size > final.size){
        # Try again
        tries   <- tries + 1
    } else {
        return(traj.mat)
    }

}

}

.sampleMolecules<-function(this, size.traj){
    final.sizes<-as.numeric(size.traj[,ncol(size.traj)])
    probs<-final.sizes/sum(final.sizes)

    while(TRUE) {
        subsamples<-as.numeric(rmultinom(
            n=1,
            size=this@sample.size,
            prob=probs 
        ) )
   
        # Multinomial sampling can result in too large subsamples:
        if(all(subsamples <= final.sizes)){
            return(subsamples)
        }
    
    }

}

