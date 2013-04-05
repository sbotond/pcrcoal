#
# Copyright (C) 2013 EMBL - European Bioinformatics Institute
#
# This program is free software: you can redistribute it
# and/or modify it under the terms of the GNU General
# Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
# Neither the institution name nor the name pcrcoal
# can be used to endorse or promote products derived from
# this software without prior written permission. For
# written permission, please contact <sbotond@ebi.ac.uk>.

# Products derived from this software may not be called
# pcrcoal nor may pcrcoal appear in their
# names without prior written permission of the developers.
# You should have received a copy of the GNU General Public
# License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

setClass(
    'PCRcoal',
    representation(
        initial.size    ="numeric",
        nr.cycles       ="numeric",
        efficiencies    ="numeric",
        sample.size     ="numeric",
        max.tries       ="numeric"
    ),
    prototype(
    ),
    sealed=TRUE
);

setGeneric("PCRcoal", function(object, ... ) standardGeneric("PCRcoal"));
setMethod(
    "PCRcoal",
    c("missing"),
    function(
        initial.size,
        nr.cycles,
        efficiencies,
        sample.size,
        max.tries
    ) {
        this<-new("PCRcoal");

        if(!missing(initial.size)){
            this@initial.size   <-initial.size;
        }

        if(!missing(nr.cycles)){
            this@nr.cycles      <-nr.cycles;
        }
        
        if(!missing(sample.size)){
            this@sample.size    <-sample.size;
        }

        if(!missing(efficiencies)){
            this@efficiencies   <-efficiencies;
        } else if (!missing(nr.cycles)){
            # The default is an efficiency of 1.0 in all cycles:
            this@efficiencies<-rep(1.0, nr.cycles)
        }

        if(!missing(max.tries)){
            this@max.tries<-max.tries
        } else {
            this@max.tries <- 100
        }

        return(this);
    }
);

setGeneric("sample.tree", function(object, ... ) standardGeneric("sample.tree"));
setMethod(
    "sample.tree",
    c("PCRcoal"),
    function(
        object
    ) {

        # Sample size trajectory:
        size.trajectory    <- .sampleTrajectory(object);

        # Sample subsample sizes:
        subsamples         <- .sampleMolecules(object, size.trajectory);
        
        # Sample tree:
        tree               <- .sampleTree(object, size.trajectory, subsamples);

        # Bless & reorder tree:
        class(tree)         <- "phylo"
        tree                <- reorder(tree)

        return(tree);
    }
);

setGeneric("sample.tnt", function(object, ... ) standardGeneric("sample.tnt"));
setMethod(
    "sample.tnt",
    c("PCRcoal"),
    function(
        object
    ) {

        # Sample size trajectory:
        size.trajectory    <- .sampleTrajectory(object);

        # Sample subsample sizes:
        subsamples         <- .sampleMolecules(object, size.trajectory);
        
        # Sample tree:
        tree               <- .sampleTree(object, size.trajectory, subsamples);

        # Bless & reorder tree:
        class(tree)         <- "phylo"
        tree                <- reorder(tree)

        return(
            list(
                "phylo"           = tree,
                "trajectories"    = size.trajectory,
                "subsamples"      = subsamples
            )
        );
    }
);

setGeneric("sample.trs", function(object, ... ) standardGeneric("sample.trs"));
setMethod(
    "sample.trs",
    c("PCRcoal"),
    function(
        object
    ) {

        # Sample size trajectory:
        size.trajectory    <- .sampleTrajectory(object);

        # Sample subsample sizes:
        subsamples         <- .sampleMolecules(object, size.trajectory);
        
        return(
            list(
                "trajectories"    = size.trajectory,
                "subsamples"      = subsamples
            )
        );
    }
);


.checkInput<-function(object){
   

    # Missing params:
    if(length(object@initial.size) == 0){
        stop("Cannot sample trajectory: initial size is missing!")
    }
    else {
        .check.integer(object@initial.size,"initial size")
    }

    if(length(object@nr.cycles) == 0){
        stop("Cannot sample trajectory: the number of cycles is missing!")
    }
    else {
        .check.integer(object@nr.cycles,"the number of cycles!")
    }

    if(length(object@sample.size) == 0){
        stop("Cannot sample trajectory: the sample size is missing!")
    }
    else {
        .check.integer(object@sample.size,"the sample size")
    }

    if(length(object@efficiencies) == 0){
        stop("Cannot sample trajectory: the efficiencies are missing!")
    }
    else if (!is.numeric(object@efficiencies)){
        stop("Cannot sample trajectory: the efficiencies must be numeric!")
    }
    else if (any(object@efficiencies > 1.0)){
        stop("Cannot sample trajectory: the efficiencies must be <= 1!")
    }

    .check.integer(object@max.tries,"max.tries") 
    
    # Dangerous situation by experience, causes integer overflow:
    if(all(object@efficiencies == 1.0) && (object@nr.cycles > 31)){
        stop("Simulation failed: all efficiencies are 1.0 and the number of cycles is larger than 31!");
    }
}

.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
    abs(x -round(x)) < tol
}

.check.integer <- function (obj, name){
    if(!.is.wholenumber(obj)){
       stop(paste("\nCannot sample trajectory: the value of ",name," slot (",obj,") is not an integer!\n\n", sep=""))
    }
}


