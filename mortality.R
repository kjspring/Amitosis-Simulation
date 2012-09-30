
#x1 <- matrix(data=rep(100, 6),ncol=1)

ip1 <-  read.csv("~/Dropbox/Research/data/bootstrap/REML/p1i.csv")
ip2 <- read.csv("~/Dropbox/Research/data/bootstrap/REML/p2i.csv")
ief <- read.csv("~/Dropbox/Research/data/bootstrap/REML/EFi.csv")
i18s <- read.csv("~/Dropbox/Research/data/bootstrap/REML/18Si.csv")

############# In this version, I added to convert 'lineages' with 0 value for copy number
############# to NA and then omitting the NA in the variance calculations
############# V5.4, this saves the mean instead of only the 1 sample
############# V5.5, added stabalizing selection to the simulation
############# V5.5.2, this allows the random seleciton of input data from the experimental
############# V5.5.3, this changes the generation time to 281 generations from initial to final
#############           gen / hour is 0.0934 generations per hour
############# V5.5.4, makes it so if a cell line is lost (selects a cell with 0 CN) it will replace it with the last one.
############# V5.5.5, runs 12 lines in parallel but seeds the 12 lines from 6
############# V6.0, measures how many of the 1000 replicates become completely extinct
############# V6.0.1 this takes 60 random numbers from a normal distribution from the original data to sim cell death

            regen <- function(x) {
                a <- which(x[1:6,1] <= 0)
                x[a] <- x[a+6]
                return(x)
                }

    sim_repo = function( x1, G=12, k=1, omega, theta=theta ) {

            # x1 is the list of copy numbers for a somatic chromosome
            # G is the number of generations, default is 12
            # k is the transfer size, default is 1

            xn <- x1
           
            # this loop does the replication for each cell in each generation
            for ( pop in 1:G ) { 
            
                # number of generations.  This is a count for the for loop
                pop <- 1 + pop 
                
                # double the somatic chromosomes for replication
                dup <- xn * 2 
                z <- matrix(rbinom(n=rep(1,length(dup)),size = as.vector(dup),prob = 0.5),nrow = nrow(dup)) 
                # amount of somatic c hromosomes distributed to one of the daughter cells
                z1 <- dup - z 
                
                # the other daughter cells receives the remainder somatic chromosomes
                xn <- cbind(z, z1) # put both in a matrix
                
                # Selection
                W <- exp( -(1/2)*( ( ( xn - theta ) / theta ) ^2 / omega ) )
                Z <- matrix(rbinom(nrow(W) * ncol(W), 1, W), nrow=nrow(W), ncol=ncol(W) ) 
                xn <- ifelse ( Z == 0, 0, xn )                
                }

            # replace 0 with NA
            xn[which(xn==0)] <- NA

            # replace any rows with only NA with 0
            xn[rowSums(is.na(xn))==ncol(xn), ] <- 0

            # find the mean of each row.  Simulates using real-time PCR
            x.mean <- round(matrix(apply(xn, 1, mean, na.rm=TRUE) ) )
            x.mean <- matrix(regen(x.mean)[1:6,1],ncol=1)

      
            # randomly select one living cell from the sample and puts it all in 1 column
            #xn <- matrix(apply(xn, 1, function(x){sample(x[!is.na(x)], size = k)}), ncol = k) # this has a bug, if there is only one value it takes a sample from i:x and not the individual
                
            xn <- matrix(apply(xn, 1, function(x)
                {if (length(x[!is.na(x)]) > 1) 
                { sample(x[!is.na(x)], size = k) } 
                else x[!is.na(x)] }), ncol=k)

            # if there is a zero result, use its pair
            #if (any(xn[1:6,1]) == 0) {
                xn <- regen(xn)
             #   }

            #output a list of the random selection and the mean of each cell line
            xn <- list(xn,x.mean)
            xn
            
    }

    # The function outputs the difference between the intial variance component between
    # 'cell lines' with the final variance after t number of transfers
    
    sim_exp = function( gene, G=12, k=1, t=23, h=1000, omega ) {

        # Add an if else statement that if rnorm selects a 0 then redo

            if (gene == "P1") {
                  #ip1 <-  read.csv("~/Dropbox/Research/data/bootstrap/REML/p1i.csv")
                  ip1mean <- mean(ip1$CNPC)
                  ip1sd <-sd(ip1$CNPC)
                  x1 <- as.matrix(round(sample(rnorm(1000,ip1mean,ip1sd),6)))
                  x1 <- matrix(rep(x1, 2), ncol=1)
                  theta <- rowMeans(x1)
                  while (any(is.na(x1))) {
                    ip2mean <- mean(ip1$CNPC)
                    ip2sd <- sd(ip1$CNPC)
                    x1 <- as.matrix(round(sample(rnorm(1000,ip1mean,ip1sd),6)))
                    x1 <- matrix(rep(x1, 2), ncol=1)
                    x1[which(x1<=0)] <- NA
                    theta <- rowMeans(x1)
                    }
                  }
                  
            if (gene == "P2") {
                  #ip2 = read.csv("~/Dropbox/Research/data/bootstrap/REML/p2i.csv")
                  ip2mean <- mean(ip2$CNPC)
                  ip2sd <- sd(ip2$CNPC)
                  x1 <- as.matrix(round(sample(rnorm(1000,ip2mean,ip2sd),60)))
                  x1 <- matrix(rep(x1, 2), ncol=1)
                  x1[which(x1<=0)] <- NA
                  theta <- rowMeans(x1)
                  while (any(is.na(x1))) {
                    ip2mean <- mean(ip2$CNPC)
                    ip2sd <- sd(ip2$CNPC)
                    x1 <- as.matrix(round(sample(rnorm(1000,ip2mean,ip2sd),60)))
                    x1 <- matrix(rep(x1, 2), ncol=1)
                    x1[which(x1<=0)] <- NA
                    theta <- rowMeans(x1)
                    }
                  }
                  
            if (gene == "EF") {
                  #ief <- read.csv("~/Dropbox/Research/data/bootstrap/REML/EFi.csv")
                  iefmean <- mean(ief$CNPC)
                  iefsd <- sd(ief$CNPC)
                  x1 <- as.matrix(round(sample(rnorm(1000,iefmean,iefsd),6)))
                  x1 <- matrix(rep(x1, 2), ncol=1)
                  theta <- rowMeans(x1)
                  while (any(is.na(x1))) {
                    ip2mean <- mean(ief$CNPC)
                    ip2sd <- sd(ief$CNPC)
                    x1 <- as.matrix(round(sample(rnorm(1000,iefmean,iefsd),6)))
                    x1 <- matrix(rep(x1, 2), ncol=1)
                    x1[which(x1<=0)] <- NA
                    theta <- rowMeans(x1)
                    }
                  }

            if (gene == "RDNA18S") {
                  #i18s <- read.csv("~/Dropbox/Research/data/bootstrap/REML/18Si.csv")
                  i18Smean <- mean(i18s$CNPC)
                  i18Ssd <- sd(i18s$CNPC)
                  x1 <- as.matrix(round(sample(rnorm(1000,i18Smean,i18Ssd),6)))
                  x1 <- matrix(rep(x1, 2), ncol=1)
                  theta <- rowMeans(x1)
                  while (any(is.na(x1))) {
                    ip2mean <- mean(i18s$CNPC)
                    ip2sd <- sd(i18s$CNPC)
                    x1 <- as.matrix(round(sample(rnorm(1000,ip2mean,i18ssd),6)))
                    x1 <- matrix(rep(x1, 2), ncol=1)
                    x1[which(x1<=0)] <- NA
                    theta <- rowMeans(x1)
                    }
                  }
                  
        x <- x1
        
        # preallocate the xn matrix
        xn <- matrix(NA, nrow(x1), t) # create a matrix

        
        # preallocate the xn.mean matrix
        xn.mean <- matrix(NA, nrow(x1)/2, t)

        # preallocate the extinct vector
        #xn.extinct <- matrix(NA, nrow=h)
        
        # enter the initial data into the preallocated matrices
        xn[,1] <- x1 
        xn.mean[,1] <- x1[1:6,1]
        
        for ( l in 2:t ) {
            # create the next population transfer
            x <- sim_repo( x1=x, G=G, k=k, omega=omega, theta=theta )
            #save them in the preallocated matrix
            xn[ ,l] <- x[[1]] 
            xn.mean[,l] <- x[[2]]
            # save the next input after the transfer
            x <- as.matrix(xn[,l])
        }

        #set any 0 values in the xn.mean matrix to NA so not to calculate the variance of missing values
        xn.mean[which(xn.mean==0)] <- NA

        # How many of the lines in the last column are extinct?
        #xn.extinct <- length(which(is.na(xn.mean[,t])))
        
        #output the variance between the first and last measurements
        colvar <- matrix(apply(xn.mean,2,var, na.rm=TRUE),ncol=ncol(xn.mean))
        ivar <- colvar[,1]
        fvar <- colvar[,ncol(xn.mean)]
        deltavar <- fvar - ivar
        y <- list(deltavar, xn.mean)
        return(y)
    }            

    # This repeats the simulation 1000 times and outputs the delta variance.

    sim_1000 <- function( gene, G=12, k=1, t=23, h=1000, omega ) {
        xn <- vector(length=h)
        xn.ext <- matrix(NA, nrow=h)

        for ( i in 1:h ) {
            x <- sim_exp( gene, G, k, t, h, omega )
            xn.ext[i,] <- length(which(is.na(x[[2]][,t])))
            xn[i] <- x[[1]]
        }
        
        xn <- list(xn, xn.ext)
        return(xn)
    } 
