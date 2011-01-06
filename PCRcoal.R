##  
## Copyright 2010- Botond Sipos  
## See the package description for licensing information.   
## 

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

setGeneric("sample.tree", function(this, ... ) standardGeneric("sample.tree"));
setMethod(
    "sample.tree",
    c("PCRcoal"),
    function(
        this
    ) {

        # Sample size trajectory:
        size.trajectory    <- .sampleTrajectory(this);

        # Sample subsample sizes:
        subsamples         <- .sampleMolecules(this, size.trajectory);
        
        # Sample tree:
        tree               <- .sampleTree(this,size.trajectory, subsamples);

        # Bless & reorder tree:
        class(tree)         <- "phylo"
        tree                <- reorder(tree)

        return(tree);
    }
);

setGeneric("sample.tnt", function(this, ... ) standardGeneric("sample.tnt"));
setMethod(
    "sample.tnt",
    c("PCRcoal"),
    function(
        this
    ) {

        # Sample size trajectory:
        size.trajectory    <- .sampleTrajectory(this);

        # Sample subsample sizes:
        subsamples         <- .sampleMolecules(this, size.trajectory);
        
        # Sample tree:
        tree               <- .sampleTree(this,size.trajectory, subsamples);

        # Bless & reorder tree:
        class(tree)         <- "phylo"
        tree                <- reorder(tree)

        return(
            list(
                "phylo"           = tree,
                "trajectories"    = size.trajectory
            )
        );
    }
);

.checkInput<-function(this){
   

    # Missing params:
    if(length(this@initial.size) == 0){
        stop("Cannot sample trajectory: initial size is missing!")
    }
    else {
        .check.integer(this@initial.size,"initial size")
    }

    if(length(this@nr.cycles) == 0){
        stop("Cannot sample trajectory: the number of cycles is missing!")
    }
    else {
        .check.integer(this@nr.cycles,"the number of cycles!")
    }

    if(length(this@sample.size) == 0){
        stop("Cannot sample trajectory: the sample size is missing!")
    }
    else {
        .check.integer(this@sample.size,"the sample size")
    }

    if(length(this@efficiencies) == 0){
        stop("Cannot sample trajectory: the efficiencies are missing!")
    }
    else if(length(this@efficiencies) != this@nr.cycles){
        stop("Cannot sample trajectory: the length of the efficiency vector must be the same as the number of cycles!")
    }
    else if (!is.numeric(this@efficiencies)){
        stop("Cannot sample trajectory: the efficiencies must be numeric!")
    }
    else if (any(this@efficiencies > 1.0)){
        stop("Cannot sample trajectory: the efficiencies must be <= 1!")
    }

    .check.integer(this@max.tries,"max.tries") 
    
    # Dangerous situation by experience:
    if(all(this@efficiencies == 1.0) && (this@nr.cycles > 31)){
        stop("Simulation failed: all efficeincies are 1.0 and the number of cycles is larger than 31!");
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

