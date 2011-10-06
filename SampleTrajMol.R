#  Copyright (C) 2011 by Botond Sipos, European Bioinformatics Institute
#  sbotond@ebi.ac.uk
#
#  This file is part of the pcrcoal software for coalescent simulations of PCR reactions.
#
#  pcrcoal is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  pcrcoal is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with pcrcoal. If not, see <http://www.gnu.org/licenses/>.

# Sample a size trajectories for each sub-population 
# (descendants of a single starting molecule)
.sampleTrajectory<-function(this){

# Check if everything is sane about the object:
.checkInput(this)

# Recycle efficiency vector:
eff <- rep(this@efficiencies, length.out=this@nr.cycles)

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
            traj.mat[i,j]<- prev.size + rbinom(n=1, size=prev.size,prob=eff[j-1]) 
        
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

