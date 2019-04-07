surface_comp_index <-
function (X.trees,Y.trees,X.Plot,Y.Plot,R.Comp) {

  n.tree <- length(X.trees)
  S <- rep(NA,n.tree)
  for (i in 1:n.tree) {
    x.tree <- X.trees[i]
    y.tree <- Y.trees[i]
  
    ####################
    # Tree in the kernel
    if (x.tree >= R.Comp & x.tree <= (X.Plot-R.Comp) & y.tree >= R.Comp & y.tree <= (Y.Plot-R.Comp)) {
      S[i] <- pi*R.Comp^2
    }
    else {
      #################
      # Tree in the top
      if (x.tree >= R.Comp & x.tree <= (X.Plot-R.Comp) & y.tree > (Y.Plot-R.Comp)) {
        # alpha coordinates
        y.alpha <- Y.Plot
        x.alpha <- x.tree-sqrt(R.Comp^2-(y.alpha-y.tree)^2)
        # m coordinates
        x.m <- x.tree
        y.m <- Y.Plot
        # length of the arrow
        a <- sqrt((y.tree-y.m)^2+(x.tree-x.m)^2)
        # length of the demi-rope
        d <- sqrt((x.alpha-x.m)^2)
        # theta angle
        theta <- 2*asin(d/R.Comp)
        # surface of the slice
        S.s <- pi*R.Comp^2*(theta/(2*pi))
        # surface of the triangle
        S.t <- a*d
        # final surface
        S[i] <- (pi*R.Comp^2)-(S.s-S.t)
      }
      ####################
      # Tree in the bottom
      if (x.tree >= R.Comp & x.tree <= (X.Plot-R.Comp) & y.tree < R.Comp) {
        # alpha coordinates
        y.alpha <- 0
        x.alpha <- x.tree-sqrt(R.Comp^2-(y.alpha-y.tree)^2)
        # m coordinates
        x.m <- x.tree
        y.m <- 0
        # length of the arrow
        a <- sqrt((y.tree-y.m)^2+(x.tree-x.m)^2)
        # length of the demi-rope
        d <- sqrt((x.alpha-x.m)^2)
        # theta angle
        theta <- 2*asin(d/R.Comp)
        # surface of the slice
        S.s <- pi*R.Comp^2*(theta/(2*pi))
        # surface of the triangle
        S.t <- a*d
        # final surface
        S[i] <- (pi*R.Comp^2)-(S.s-S.t)
      }
      ##################
      # Tree in the left
      if (x.tree < R.Comp & y.tree >= R.Comp & y.tree <= (Y.Plot-R.Comp)) {
        # alpha coordinates
        x.alpha <- 0
        y.alpha <- y.tree-sqrt(R.Comp^2-(x.alpha-x.tree)^2)
        # m coordinates
        x.m <- 0
        y.m <- y.tree
        # length of the arrow
        a <- sqrt((y.tree-y.m)^2+(x.tree-x.m)^2)
        # length of the demi-rope
        d <- sqrt((y.alpha-y.m)^2)
        # theta angle
        theta <- 2*asin(d/R.Comp)
        # surface of the slice
        S.s <- pi*R.Comp^2*(theta/(2*pi))
        # surface of the triangle
        S.t <- a*d
        # final surface
        S[i] <- (pi*R.Comp^2)-(S.s-S.t)
      }
      ###################
      # Tree in the right
      if (x.tree > (X.Plot-R.Comp) & y.tree >= R.Comp & y.tree <= (Y.Plot-R.Comp)) {
        # alpha coordinates
        x.alpha <- X.Plot
        y.alpha <- y.tree-sqrt(R.Comp^2-(x.alpha-x.tree)^2)
        # m coordinates
        x.m <- X.Plot
        y.m <- y.tree
        # length of the arrow
        a <- sqrt((y.tree-y.m)^2+(x.tree-x.m)^2)
        # length of the demi-rope
        d <- sqrt((y.alpha-y.m)^2)
        # theta angle
        theta <- 2*asin(d/R.Comp)
        # surface of the slice
        S.s <- pi*R.Comp^2*(theta/(2*pi))
        # surface of the triangle
        S.t <- a*d
        # final surface
        S[i] <- (pi*R.Comp^2)-(S.s-S.t)
      }
      #########################
      # Tree in top left corner
      if (x.tree < R.Comp & y.tree > (Y.Plot-R.Comp)) {
        if (sqrt((x.tree-0)^2+(y.tree-Y.Plot)^2) <= R.Comp) { # Circle overlaps the corner
          # alpha coordinates
          x.alpha <- 0
          y.alpha <- y.tree-sqrt(R.Comp^2-(x.alpha-x.tree)^2)
          # beta coordinates
          y.beta <- Y.Plot
          x.beta <- x.tree+sqrt(R.Comp^2-(y.beta-y.tree)^2)
          # m coordinates
          x.m <- (x.alpha+x.beta)/2
          y.m <- (y.alpha+y.beta)/2
          # length of the arrow
          a <- sqrt((y.tree-y.m)^2+(x.tree-x.m)^2)
          # length of the demi-rope
          d <- sqrt((y.alpha-y.beta)^2+(x.alpha-x.beta)^2)/2
          # theta angle
          theta <- 2*asin(d/R.Comp)
          # surface of the slice
          S.s <- pi*R.Comp^2*(theta/(2*pi))
          # surface of the first triangle
          S.ft <- a*d
          # surface of the second triangle
          S.st <- x.beta*(Y.Plot-y.alpha)/2
          # final surface
          S[i] <- S.s-S.ft+S.st
        }
        else { # Circle does not overlap the corner
          # alpha coordinates
          x.alpha <- 0
          y.alpha <- y.tree-sqrt(R.Comp^2-(x.alpha-x.tree)^2)
          # beta coordinates
          y.beta <- Y.Plot
          x.beta <- x.tree+sqrt(R.Comp^2-(y.beta-y.tree)^2)
          # length of the arrow 1
          a1 <- sqrt((x.tree-0)^2)
          # length of the demi-rope 1
          d1 <- sqrt((y.alpha-y.tree)^2)
          # theta 1 angle
          theta1 <- 2*asin(d1/R.Comp)
          # surface of the slice 1
          S.s1 <- pi*R.Comp^2*(theta1/(2*pi))
          # surface of the triangle 1
          S.t1 <- a1*d1
          # length of the arrow 2
          a2 <- sqrt((y.tree-Y.Plot)^2)
          # length of the demi-rope 2
          d2 <- sqrt((x.beta-x.tree)^2)
          # theta 2 angle
          theta2 <- 2*asin(d2/R.Comp)
          # surface of the slice 2
          S.s2 <- pi*R.Comp^2*(theta2/(2*pi))
          # surface of the triangle 2
          S.t2 <- a2*d2
          # final surface
          S[i] <- (pi*R.Comp^2)-(S.s1-S.t1)-(S.s2-S.t2)
        }
      }
      ##########################
      # Tree in top right corner
      if (x.tree > (X.Plot-R.Comp) & y.tree > (Y.Plot-R.Comp)) {
        if (sqrt((x.tree-X.Plot)^2+(y.tree-Y.Plot)^2) <= R.Comp) { # Circle overlaps the corner	
          # alpha coordinates
          x.alpha <- X.Plot
          y.alpha <- y.tree-sqrt(R.Comp^2-(x.alpha-x.tree)^2)
          # beta coordinates
          y.beta <- Y.Plot
          x.beta <- x.tree-sqrt(R.Comp^2-(y.beta-y.tree)^2)
          # m coordinates
          x.m <- (x.alpha+x.beta)/2
          y.m <- (y.alpha+y.beta)/2
          # length of the arrow
          a <- sqrt((y.tree-y.m)^2+(x.tree-x.m)^2)
          # length of the demi-rope
          d <- sqrt((y.alpha-y.beta)^2+(x.alpha-x.beta)^2)/2
          # theta angle
          theta <- 2*asin(d/R.Comp)
          # surface of the slice
          S.s <- pi*R.Comp^2*(theta/(2*pi))
          # surface of the first triangle
          S.ft <- a*d
          # surface of the second triangle
          S.st <- (X.Plot-x.beta)*(Y.Plot-y.alpha)/2
          # final surface
          S[i] <- S.s-S.ft+S.st
        }
        else { # Circle does not overlap the corner
          # alpha coordinates
          x.alpha <- X.Plot
          y.alpha <- y.tree-sqrt(R.Comp^2-(x.alpha-x.tree)^2)
          # beta coordinates
          y.beta <- Y.Plot
          x.beta <- x.tree-sqrt(R.Comp^2-(y.beta-y.tree)^2)
          # length of the arrow 1
          a1 <- sqrt((x.tree-X.Plot)^2)
          # length of the demi-rope 1
          d1 <- sqrt((y.alpha-y.tree)^2)
          # theta 1 angle
          theta1 <- 2*asin(d1/R.Comp)
          # surface of the slice 1
          S.s1 <- pi*R.Comp^2*(theta1/(2*pi))
          # surface of the triangle 1
          S.t1 <- a1*d1
          # length of the arrow 2
          a2 <- sqrt((y.tree-Y.Plot)^2)
          # length of the demi-rope 2
          d2 <- sqrt((x.beta-x.tree)^2)
          # theta 2 angle
          theta2 <- 2*asin(d2/R.Comp)
          # surface of the slice 2
          S.s2 <- pi*R.Comp^2*(theta2/(2*pi))
          # surface of the triangle 2
          S.t2 <- a2*d2
          # final surface
          S[i] <- (pi*R.Comp^2)-(S.s1-S.t1)-(S.s2-S.t2)
        }
      }
      ############################
      # Tree in bottom left corner
      if (x.tree < R.Comp & y.tree < R.Comp) {
        if (sqrt((x.tree-0)^2+(y.tree-0)^2) <= R.Comp) { # Circle overlaps the corner
          # alpha coordinates
          x.alpha <- 0
          y.alpha <- y.tree+sqrt(R.Comp^2-(x.alpha-x.tree)^2)
          # beta coordinates
          y.beta <- 0
          x.beta <- x.tree+sqrt(R.Comp^2-(y.beta-y.tree)^2)
          # m coordinates
          x.m <- (x.alpha+x.beta)/2
          y.m <- (y.alpha+y.beta)/2
          # length of the arrow
          a <- sqrt((y.tree-y.m)^2+(x.tree-x.m)^2)
          # length of the demi-rope
          d <- sqrt((y.alpha-y.beta)^2+(x.alpha-x.beta)^2)/2
          # theta angle
          theta <- 2*asin(d/R.Comp)
          # surface of the slice
          S.s <- pi*R.Comp^2*(theta/(2*pi))
          # surface of the first triangle
          S.ft <- a*d
          # surface of the second triangle
          S.st <- x.beta*y.alpha/2
          # final surface
          S[i] <- S.s-S.ft+S.st
        }
        else { # Circle does not overlap the corner
          # alpha coordinates
          x.alpha <- 0
          y.alpha <- y.tree+sqrt(R.Comp^2-(x.alpha-x.tree)^2)
          # beta coordinates
          y.beta <- 0
          x.beta <- x.tree+sqrt(R.Comp^2-(y.beta-y.tree)^2)
          # length of the arrow 1
          a1 <- sqrt((x.tree-0)^2)
          # length of the demi-rope 1
          d1 <- sqrt((y.alpha-y.tree)^2)
          # theta 1 angle
          theta1 <- 2*asin(d1/R.Comp)
          # surface of the slice 1
          S.s1 <- pi*R.Comp^2*(theta1/(2*pi))
          # surface of the triangle 1
          S.t1 <- a1*d1
          # length of the arrow 2
          a2 <- sqrt((y.tree-0)^2)
          # length of the demi-rope 2
          d2 <- sqrt((x.beta-x.tree)^2)
          # theta 2 angle
          theta2 <- 2*asin(d2/R.Comp)
          # surface of the slice 2
          S.s2 <- pi*R.Comp^2*(theta2/(2*pi))
          # surface of the triangle 2
          S.t2 <- a2*d2
          # final surface
          S[i] <- (pi*R.Comp^2)-(S.s1-S.t1)-(S.s2-S.t2)
        }
      }
      #############################
      # Tree in bottom right corner
      if (x.tree > (X.Plot-R.Comp) & y.tree < R.Comp) {
        if (sqrt((x.tree-X.Plot)^2+(y.tree-0)^2) <= R.Comp) { # Circle overlaps the corner
          # alpha coordinates
          x.alpha <- X.Plot
          y.alpha <- y.tree+sqrt(R.Comp^2-(x.alpha-x.tree)^2)
          # beta coordinates
          y.beta <- 0
          x.beta <- x.tree-sqrt(R.Comp^2-(y.beta-y.tree)^2)
          # m coordinates
          x.m <- (x.alpha+x.beta)/2
          y.m <- (y.alpha+y.beta)/2
          # length of the arrow
          a <- sqrt((y.tree-y.m)^2+(x.tree-x.m)^2)
          # length of the demi-rope
          d <- sqrt((y.alpha-y.beta)^2+(x.alpha-x.beta)^2)/2
          # theta angle
          theta <- 2*asin(d/R.Comp)
          # surface of the slice
          S.s <- pi*R.Comp^2*(theta/(2*pi))
          # surface of the first triangle
          S.ft <- a*d
          # surface of the second triangle
          S.st <- (X.Plot-x.beta)*y.alpha/2
          # final surface
          S[i] <- S.s-S.ft+S.st
        }
        else { # Circle does not overlap the corner
          # alpha coordinates
          x.alpha <- X.Plot
          y.alpha <- y.tree+sqrt(R.Comp^2-(x.alpha-x.tree)^2)
          # beta coordinates
          y.beta <- 0
          x.beta <- x.tree-sqrt(R.Comp^2-(y.beta-y.tree)^2)
          # length of the arrow 1
          a1 <- sqrt((x.tree-X.Plot)^2)
          # length of the demi-rope 1
          d1 <- sqrt((y.alpha-y.tree)^2)
          # theta 1 angle
          theta1 <- 2*asin(d1/R.Comp)
          # surface of the slice 1
          S.s1 <- pi*R.Comp^2*(theta1/(2*pi))
          # surface of the triangle 1
          S.t1 <- a1*d1
          # length of the arrow 2
          a2 <- sqrt((y.tree-0)^2)
          # length of the demi-rope 2
          d2 <- sqrt((x.beta-x.tree)^2)
          # theta 2 angle
          theta2 <- 2*asin(d2/R.Comp)
          # surface of the slice 2
          S.s2 <- pi*R.Comp^2*(theta2/(2*pi))
          # surface of the triangle 2
          S.t2 <- a2*d2
          # final surface
          S[i] <- (pi*R.Comp^2)-(S.s1-S.t1)-(S.s2-S.t2)
        }
      }
    }
  }
  return(S)
}

