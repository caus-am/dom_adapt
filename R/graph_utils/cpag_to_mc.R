cpag_to_mc <- function(P){
  # convert cPAG P into causal matrix Mc, including indirect non/causal info
  #
  # Input convention:
  # The edge marks are encoded by numbers: 0 = no edge, 1 = circle, 2 = arrowhead, 3 = tail.
  #
  # Internal convention:
  # The edge marks are encoded by numbers: 0 = no edge, 1 = tail, 2 = arrowhead, 3 = circle.
  
  n <- nrow(P);
  
  input_circle <- 1;
  input_tail <- 3;
  
  noedge <- 0;
  circle <- 3;
  arrowhead <- 2;
  tail <- 1;
  
  P[which(P==input_circle)] = -3;
  P[which(P==input_tail)] = tail;
  P[which(P==-3)] = circle;
  
  
  # first include all implicit causal relation (= Zhang rule R10, 4.2 in LoCI)
  # note: this was missing from the implemention in the original evaluation
  # as pointed out by Sofia, see Th3.1 in 'Marginal causal consistency ..'
  # loop over all unshielded non-colliders with uncov.pot.dir.paths to y
  Mimpc = Get_Implicit_Causal(P,n);
  
  # next construct reachability graph from C:  
  # compute transcausal: remove all arrowheads (-1) and circles
  Mtc = P;
  Mtc[which(P == arrowhead | P == circle)] = 0;
  Mtc = Mtc + Mimpc;
  # causal reachability matrix now becomes
  #Mc = sign(expm(Mtc)) - eye(N);  # reachable by directed path
  Mc = sign(expm(Mtc)) - diag(n);
  
  # idem, non-reachable without going against an arrowhead
  Mtnc = P;
  Mtnc[which(P == arrowhead)] = 0;     # remove arrowhead
  Mtnc[which(P == circle)] = tail;     # circles become tails
  Mtnc = (sign(expm(Mtnc)) -  diag(n)) - array(1, c(n, n));
  # non-causal reachability matrix now becomes
  Mtnc = sign(Mtnc - t(Mc)); # subtract reverse causal (not necessary for pure PAGs)

  # causal matrix is 
  Mc = Mc + Mtnc;
  
  # Convert to standard notation in this framework:
  newC <- t(Mc)
  
  list(C=newC, G=array(0, c(n, n)), Ge=array(0, c(n, n)))
}
  

cpag_to_mc.old <- function(C){
  # convert cPAG P into causal matrix Mc, including indirect non/causal info
  #
  # Old convention: C(i,j) = {-2,-1,0,1}
  # if there is a {unknown,tail,nothing,arrow} at j
  #
  # New convention:
  # The edge marks are encoded by numbers: 0 = no edge, 1 = circle, 2 = arrowhead, 3 = tail.
  #
  #
  # MC(i,j) = 1 if j cannot be an ancestor of i (arrow)
  # MC(i,j) = -1 if j must be an ancestor of i (tail)
  # MC(i,j) = -2 if unknown
  
  n <- nrow(C);
  
  #acausalp = 1 * (C == -1) + 1 * (C == -2); % reachable without going against arrowheads
  notCauses <- 1*(C==3) + 1*(C==1)
  
  #acausalp = (expm(acausalp') == 0); % transitivity (reachable -> not sure that acausal)
  notCauses <- 1* (expm(t(notCauses)) == 0)
  #acausalp = acausalp + eye(N);
  notCauses <- notCauses + diag(n)
  
  # causalp = (C == 1) .* (C' == -1); % reachable via directed edges
  # reachable via directed edges
  causes <- (1*(C==2)) * (1*(t(C)==3))
  
  # transitivity
  # causalp = (expm(causalp) - eye(N) > 0); % transitivity
  causes <- 1*(expm(causes) - diag(n) > 0)
  
  #MC = acausalp - causalp;
  #MC = MC + -2 * (MC == 0);
  MC <- notCauses - causes
  MC <- MC + -2 * (MC==0)
  
  # MC(i,j) = 1 if j cannot be an ancestor of i (arrow)
  # MC(i,j) = -1 if j must be an ancestor of i (tail)
  # MC(i,j) = -2 if unknown
  
  # Convert to standard notation in this framework:
  newC <- MC
  newC[which(MC == -2)] <-  0
  newC[which(MC == 1)]  <- -1
  newC[which(MC == -1)] <-  1
  
  list(C=newC, G=array(0, c(n, n)), Ge=array(0, c(n, n)))
}

test_cpdag_to_mc <- function(){
  # 1 <-> 2
  C <- matrix(c(0,2,2,0), 2)
  res <- matrix(c(-1,-1,-1,-1), 2)
  if (any(cpag_to_mc(C)$C!=res)) cat("ERROR 1!")
  
  # 1 o-o 2
  C <- matrix(c(0,1,1,0), 2)
  res <- matrix(c(-1,0,0,-1), 2)
  if (any(cpag_to_mc(C)$C!=res)) cat("ERROR 2!")
  
  # 1 -> 2
  C <- matrix(c(0,2,3,0), 2)
  res <- matrix(c(-1,1,-1,-1), 2)
  if (any(cpag_to_mc(C)$C!=res)) cat("ERROR 3!")
  
  # 1 o-> 2
  C <- matrix(c(0,1,2,0), 2)
  res <- matrix(c(-1,-1,0,-1), 2)
  if (any(cpag_to_mc(C)$C!=res)) cat("ERROR 4!")
  

  C = c(0,1,0,0,2,1, 1,0,1,0,0,0, 0,2,0,3,0,0, 0,0,2,0,0,2, 3,0,0,0,0,1, 1,0,0,3,0,0);
  C = matrix(C,nrow = 6,ncol = 6)
  C = t(C)

  res = c(-1, 0, 0, 1, 0, 0,
          0, -1, 0, 1,-1, 0,
         -1, -1,-1, 1,-1,-1,
         -1, -1,-1,-1,-1,-1,
          1,  0, 0, 1,-1, 0,
          0,  0, 0, 1,-1,-1);
  res = matrix(res,nrow = 6,ncol = 6)
  if (any(cpag_to_mc(C)$C!=res)) cat("ERROR 5!")
  
  C = c(0,1,0,0,2,1, 1,0,1,0,1,0, 0,2,0,3,0,0, 0,0,2,0,0,2, 3,1,0,0,0,1, 1,0,0,3,1,0);
  C = matrix(C,nrow = 6,ncol = 6)
  C = t(C)
  
  res = c(-1, 0, 0, 1,-1, 0,
          0, -1, 0, 1,0, 0,
          -1, -1,-1, 1,-1,-1,
          -1, -1,-1,-1,-1,-1,
          1,  0, 0, 1,-1, 0,
          0,  0, 0, 1,0,-1);
  
  res = matrix(res,nrow = 6,ncol = 6)
  if (any(cpag_to_mc(C)$C!=res)) cat("ERROR 6!")
}
      


Get_Implicit_Causal <- function (P,n) {
  
  noedge <- 0;
  circle <- 3;
  arrowhead <- 2;
  tail <- 1;

  Mimpc = array(0, c(n,n));
  # = orientation rule R10 from FCI: double uncovered pdp on unsh.noncollider
  # if x (noedge) y, q --> y <-- z, u1=(x,(s,..,)z) and u2=(x,(q,..,)r) u.p.d. paths 
  # then if s (possibly z) and q (possibly r) not adjacent, then orient x --> y
  # NOTE: the final arcs z --> y <-- v necessarily result from R9
  
  for (y in 1:n){
    # get all nodes with z --> y 
    # relax to also allow z o-> y from non-pure PAGs
    # Z = find( P(y,:) == 2 & (P(:,y)' == 1 | P(:,y)' == 3));
    Z = which( P[y,] == arrowhead & (t(P[,y]) == tail | t(P[,y]) == circle));
    
    nZ = length(Z);
    if (nZ < 2){ next; }
    # get all other x with no edge to y 
    # NOTE: this is the only difference with rule R10: x o-> y
    #X = mysetdiff(find( P(y,:) == 0 & P(:,y)' == 0),y);
    X = which(P[y,] == noedge & t(P[,y]) == noedge);
    X = X[!X==y];
    
    for (x in t(X)){
      # get all nodes with tail/circle adjacent to x
      # note: no arrowhead, as e.g. s *-> x would imply x --> p,
      # and so on towards y, and then back to x via R8a, and so x --> y (by transitivity)
      
      S = which( (P[x,] == tail | P[x,] == circle) & (t(P[,x]) == arrowhead | t(P[,x]) == circle) );
      nS = length(S);
      if (nS < 2) { next; }
      # find first uncovered pdp, if found try second
      for (i in 1:(nS-1)){
        s = S[i];         
        # find all nodes in S not adjacent to s (used in up2)
        Q = S[which(P[s,S] == 0)];
        Q = Q[!Q==s];
        nQ = length(Q);
        if (nQ == 0){ next;}
        # ok: try first upd path x - S(i) .. Z(j)
        # handle s == z?
        for (j in 1:nZ){     # note: not nZ-1
          z = Z[j];
          # find uncovered pdp x -o s - .. z
          # note: we cannot look for just 'any' path (as in R9), as it is
          # about a specific combination of two disjoint paths 
          # store beginning of possible uncovered possibly directed paths [x (-/o)-(o/>) s (..) y]
          up1 <- list();
          up1[[1]] = c(x,s,z);

          # process edge
          while (length(up1) > 0 ) {
            u1 = up1[[1]][length(up1[[1]])-2];   # yes: handle s==z (-> go to up2
            v1 = up1[[1]][length(up1[[1]])-1];
            # get possible extension nodes: u o-> v o-> w, with w not adjacent to u)
            idxW1 = (P[v1,]  == tail | P[v1,] == circle ) & (t(P[,v1]) == arrowhead | t(P[,v1]) == circle ) & (P[u1,] == 0);
            # check if we're there: reach target z (no additional checks needed)  
            
            if ((idxW1[z] > 0) | (s == z)) {
              # yes: uncovered possibly directed path: try second path!
              # get all other nodes in Z
              R = Z[!Z==z];
              # now loop for second path
              for (k in 1:nQ) {            # note: not nQ-1
                q = Q[k];
                for (l in 1:(nZ-1)){     # nR = nZ - 1
                  r = R[l];
                  up2 <- list();
                  up2[[1]] = c(x,q,r);   # note: up2{one}

                  # process edge (identical to first path
                  
                  while(length(up2) > 0 ) {
                    u2 = up2[[1]][length(up2[[1]])-2]; # avoid overwriting {u1,v1}
                    v2 = up2[[1]][length(up2[[1]])-1];
                    # get possible extension nodes: u o-> v o-> w, with w not adjacent to u)
                    
                    idxW2 = (P[v2,]  == tail | P[v2,] == circle ) & (t(P[,v2]) == arrowhead | t(P[,v2]) == circle ) & (P[u2,] == 0);
                    
                    # check if we're there: reach target z (no additional checks needed)
                    if ((idxW2[r] > 0) || (q == r)) { 
                      # yes: second uncovered possibly directed path!
                      Mimpc[x,y] = 1;
                      break;    # go to next unchecked edge (outer while loop)
                    }
                    else {
                      # no: exclude nodes already on the path (and/or the
                      # first path? ... yes: otherwise other rule should
                      # have triggered already
                      idxW2[up2[[1]]] = 0;
                      W2 = which(idxW2 > 0);
                      # extend path with 
                      for (w2 in t(W2)){
                        up2[[length(up2)+1]] = c(up2[[1]][1:(length(up2[[1]])-1)], w2, up2[[1]][length(up2[[1]])]);
                      }
                      # try next path
                      if (length(up2) >= 2){
                        up2 = up2[2:length(up2)];
                      } else {
                        up2 = NULL;
                      }
                    }
                  } # while up2
                  if (Mimpc[x,y] == 1){ break }
                 } # for l
                 if (Mimpc[x,y] == 1){ break }
                } # for k
              } else {  
               # z not reached in idxW1
               # exclude nodes already on the path 
                idxW1[up1[[1]]] = 0;
                W1 = which(idxW1 > 0);
                # extend path with 
                for (w1 in t(W1)) {
                  up1[[length(up1)+1]] = c(up1[[1]][1:(length(up1[[1]])-1)],w1,up1[[1]][length(up1[[1]])]);
                }
                
              } #if (z in idxW1)
              # try next path
            
            if (length(up1) >= 2){
              up1 = up1[2:length(up1)];
            } else {
              up1 = NULL;
            }
            
            if (Mimpc[x,y] == 1){ break }
           } # while up1
          if (Mimpc[x,y] == 1) { break }
        } # for j
        if (Mimpc[x,y] == 1) {break }
      }  # for i 
    } # for x
  } # for y
  Mimpc
}