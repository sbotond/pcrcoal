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

.sampleTree<-function(this, size.traj, subsam){

      # Initialize pools
      tmp       <-.initPools(subsam) 
      pools     <-tmp$pools
      max.node  <-tmp$max
      rm(tmp)
      
      # Number of the leafs in the final tree:
      nr.leafs  <-max.node

      # "Declare" phylo object elements
      edge.from <-integer()
      edge.to   <-integer()
      edge.len  <-integer()

      # Base node pool:
      base.pool <-list()
      tp        <-list()
           
      for(s in names(pools) ){

           nr   <- as.numeric(s)
           tp   <- .traceSubsam(pools[[s]], as.numeric(size.traj[nr,]), max.node)

           # Update maximum node id:
           max.node <- tp$max.node

           # Update tree elements:
           edge.from    <-  c(edge.from, tp$tree$edge.from)  
           edge.to      <-  c(edge.to, tp$tree$edge.to)  
           edge.len     <-  c(edge.len, tp$tree$edge.len)  

           # Collect base node:
           base.pool[[tp$base.node]]    <- tp$base.bl
           
      }
      
      # Coalesce base nodes: 
      tp   <-  .coalBase(base.pool, max.node, subsam, this)
     
      max.node  <- tp$max.node 
      # Update tree elements:
      edge.from    <-  c(edge.from, tp$edge.from)  
      edge.to      <-  c(edge.to, tp$edge.to)  
      edge.len     <-  c(edge.len, tp$edge.len)  

      # Create the APE phylo object:
      phylo<-list()

      phylo$edge            <-cbind(edge.from, edge.to)
      phylo$edge.length     <-edge.len
      phylo$Nnode           <-(max.node - nr.leafs)
      phylo$tip.label       <-paste("m",1:nr.leafs,sep="")

      phylo                 <-.mapInternalNodes(phylo, max.node, nr.leafs)

      # Check the sanity of the branch lengths:
      if( any(phylo$edge.length > length(size.traj)) ){
        stop("\nOne of the branch lengths is higher than the number of cycles!\nThe simulation is flawed!\n\n");
      }

      return(phylo)
}

.mapInternalNodes<-function(p, max.node, nr.leafs){

    # Check for weirdness:
    if(p$Nnode < 1) { stop("\n\nSimulation resulted in a tree with a single node (which is invalid)!\n\n") }
    if(p$Nnode == 1){ return(p) }

    # APE requires the root node to have the id "max.leaf + 1"
    # , so the internal nodes must be mapped:
    int.nodes   <- (nr.leafs +1):max.node
    rev.nodes   <- rev(int.nodes)

    node.map    <- list() 
   
    # Build node map: 
    for(i in 1:length(int.nodes)) {
        node.map[[ as.character(int.nodes[i]) ]]    <- rev.nodes[i]
    }

	rm(int.nodes)
	rm(rev.nodes)

    d   <- dim(p$edge) 

    # Substitute nodes:
    p$edge<-as.numeric(p$edge)

    for(i in 1:length(p$edge) ) {
        replacement <- node.map[[ as.character(p$edge[i]) ]]
        if(!is.null(replacement)){
            p$edge[i]   <- replacement
        }
    }

    dim(p$edge) <-  d

    return(p)
}

.coalBase<-function(base.pool, max.node, subsample, this){
    res <- list()

    res$edge.from<-integer()
    res$edge.to  <-integer()
    res$edge.len <-integer()

    while(length(base.pool) > 1) {
       # Sample two nodes:
       coal.nodes<-sample(names(base.pool),2,replace=FALSE)  

       # Create parent node:
       max.node <- max.node + 1

       # Create the first edge: 
       res$edge.from    <-  c(res$edge.from, max.node) 
       res$edge.to      <-  c(res$edge.to, as.numeric(coal.nodes[1]) )
       
       # Create the second edge: 
       res$edge.from    <-  c(res$edge.from, max.node) 
       res$edge.to      <-  c(res$edge.to, as.numeric(coal.nodes[2]) )

       # Set the branch lengths:
       res$edge.len <-  c(res$edge.len, base.pool[[ coal.nodes[1] ]])
       res$edge.len <-  c(res$edge.len, base.pool[[ coal.nodes[2] ]])

       # Add the parent node to pool with replication count zero:
       base.pool[[ as.character(max.node) ]] <- 0

       # Remove child nodes:
       base.pool[[ coal.nodes[1] ]]    <- NULL
       base.pool[[ coal.nodes[2] ]]    <- NULL

    }

    # Sanity check:
    if(length(base.pool) != 1) {
        stop("Error when coalescing base nodes!");
    }

    res$max.node    <- max.node
    
    # Deal with the last node:
    last.node   <- names(base.pool)[1]

    # What if last node has branch length? 
    if(base.pool[[last.node]] > 0) {

       # This implies that all but one subsample has size zero!
       if( !any(subsample == this@sample.size) ){
        stop("The final node has replication count > 0, yet more than one subsample has non-zero size!")
       }

       # Create the "ultimate node":
       ultimate.node    <- max.node + 1
       # Create edge:
       res$edge.from    <-  c(res$edge.from, ultimate.node) 
       res$edge.to      <-  c(res$edge.to, as.numeric(last.node) )

       # Set branch lengths:
       res$edge.len   <-  c(res$edge.len, base.pool[[ last.node ]])

       res$max.node  <- ultimate.node
    }

    return(res)
}

.traceSubsam<-function(pool, size.traj, max.node){
    res<-list()
    tree<-list()
    
    # "Declare" phylo object elements
    tree$edge.from <-integer()
    tree$edge.to   <-integer()
    tree$edge.len  <-integer()
    
    # Set initial number of molecules to initial pool size:
    Ni  <- length(pool) 

    # Iterate back over size trajectory:
    for( cycle in (length(size.traj) - 1):1 ) {

       # Sample number of molecules synthetized in the
       # current cycle: 
       Ri   <- .sampleRi(cycle, size.traj, Ni)

       # Sample the number of coalescent events:
       Li   <- .sampleLi(cycle, size.traj, Ni, Ri)

       # Update the number of the nodes:
       Ni   <- Ni - Li

       tmp<-list(
        new.nodes=character()
       )

       # Coalesce nodes if Li > 0
       if (Li > 0){
            tmp<-.coalNodes(pool, Li, tree, max.node);

            # Update data structures:
            pool        <- tmp$pool
            max.node    <- tmp$max.node
            tree        <- tmp$tree
       }

       # Update synthesis count:
       Ri   <- Ri - Li

       # Distribute remaining synthesis count: 
       snodes   <-  sample(
            setdiff( names(pool), tmp$new.nodes),
            Ri,
            replace=FALSE
        ) 

        for(sn in snodes){
            pool[[sn]]   <- pool[[sn]] + 1
        }
        
        # Set synthesis count ot zero:
        Ri  <- 0

    } # for cycle

    # Check if the base node is ok:
    if(length(pool) != 1) {
        stop("More than one base node in subsample!")
    }

    res$tree        <- tree
    res$max.node    <- max.node
    res$base.node   <- names(pool)[1]   # This is not necessarily max.node
    res$base.bl     <- pool[[1]]

    return(res)
}

.coalNodes<-function(pool, Li, tree, max.node){
    res<-list()

    new.nodes   <- integer()

    # Sample nodes to coalesce:
    coal.nodes  <- sample(names(pool), (Li * 2), replace=FALSE )

    # Iterate by doublets:
    for(i in seq( from=1, to=(length(coal.nodes) - 1), by=2 ) ){
        
        # Create parent node:
        max.node <- max.node + 1
        
        # Add to new nodes list:
        new.nodes<-c(new.nodes, as.character(max.node))

        # Append the first edge:
        tree$edge.from  <-  c(tree$edge.from, max.node)
        tree$edge.to    <-  c(tree$edge.to, as.numeric(coal.nodes[i]))
        
        # Append the second edge:
        tree$edge.from  <-  c(tree$edge.from, max.node)
        tree$edge.to    <-  c(tree$edge.to, as.numeric(coal.nodes[i + 1]))

        # Set the branch lengths (note, that one of them is replicated!):
        bl.plus         <- sample(c(1, 0), 2, replace=FALSE )

        tree$edge.len   <- c(tree$edge.len, (pool[[ coal.nodes[i] ]]     + bl.plus[1]) )
        tree$edge.len   <- c(tree$edge.len, (pool[[ coal.nodes[i + 1] ]] + bl.plus[2]) )

        # Add parent node to pool with bl count of zero:
        pool[[ as.character(max.node) ]]    <- 0
        
        # Remove children nodes:
        pool[[ coal.nodes[i] ]]         <- NULL
        pool[[ coal.nodes[i + 1] ]]     <- NULL
        
    }

    res$pool        <- pool
    res$tree        <- tree
    res$max.node    <- max.node
    res$new.nodes   <- new.nodes

    return(res)
}

.sampleRi<-function(cycle, size.traj, Ni){
    # Cycle index reversed compared to the article!!
        
    sdiff<-size.traj[cycle + 1] - size.traj[cycle]
    ri<-rhyper(
        nn  = 1,       
        m   = sdiff,                # white balls (synthetized)
        n   = size.traj[cycle],     # black balls (not synthetized)
        k   = Ni                    # balls drawn (total molecules in cycle i)
    )

    if(is.nan(ri)){
        stop(paste("\nSomething is wrong! Ri is NaN in cycle: ", cycle,"\n\n"),sep="");
    }

    return(ri)
}

.sampleLi<-function(cycle, size.traj, Ni, Ri){
    li<-rhyper(
        nn  =1,
        m   = (Ni - Ri),                        # white balls (not-syntheetized but in subsample)
        n   = (size.traj[cycle] - (Ni - Ri)),   # black balls (not synthetized, not in subsample)
        k   = Ri                                # balls drawn (total molecules synthetized)
    ) 
    
    if(is.nan(li)){
        stop(paste("Something is wrong! Li is NaN in cycle:", cycle));
    }

    return(li)
}

.initPools<-function(subsam){
    pools   <-list()
    max     <-0

    # Iterate over subsamples:
    for(i in seq(along.with=subsam)){

        if(subsam[i] > 0){
            tmp<-.initOnePool(subsam[i], max)

            if(length(tmp$pool) != subsam[i]){
                stop("Faulty pool initialization!")
            }

            max<- tmp$max
            pools[[ as.character(i) ]] <- tmp$pool
        }

    }
    
    return(
        list(
            "pools" = pools,
            "max"   = max
        )
    )
}

.initOnePool<-function(size, max){
    res         <- list()
    res$pool    <- list()
    new.max     <- size + max

    for(i in (max+1):(new.max) ){
       # Initilaize replication counter to zero
       # for the leaf nodes:
       res$pool[[ as.character(i) ]] <- 0
    }

    res$max <-  new.max
    return(res)
}
