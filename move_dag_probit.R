n.edge = function(A){
  length(which(A[lower.tri(A)] == 1 | t(A)[lower.tri(A)] == 1))
}

names   = c("action","test","x","y")
actions = c("id","dd","rd")

# types are then indexed by (1,2,3,4,5)

move = function(A, q = q){

  # A: adjacency matrix of the DAG
  # q: number of vertices

  # Output: 
  # A direct successor DAG
  # the (estimated) number of direct successors
  
  A_na = A_na_na = A
  
  A_na[1,] = NA
  
  A_na_na[1,] = A_na_na[,1] = NA
  
  diag(A_na) = diag(A_na_na) = NA
  
  id_set = c()
  dd_set = c()
  rd_set = c()
  
  # set of nodes for id
  
  set_id = which(A_na == 0, TRUE)
  
  if(length(set_id) != 0){
    id_set = cbind(1, rbind(set_id, set_id[,1:2]))
  }
  
  # set of nodes for dd
  
  set_dd = which(A == 1, TRUE)
  
  if(length(set_dd != 0)){
    dd_set = cbind(2, set_dd)
  }
  
  # set of nodes for rd
  
  set_rd = which(A_na_na == 1, TRUE)
  
  if(length(set_rd != 0)){
    rd_set = cbind(3, set_rd)
  }
  
  O = rbind(id_set, dd_set, rd_set)

  repeat {
    
    i = sample(dim(O)[1],1)
    
      act_to_exe  = paste0(actions[O[i,1]],"(A=A,c(",as.vector(O[i,2]),",",as.vector(O[i,3]),"))")
      A_succ      = eval(parse(text = act_to_exe))
      act_to_eval = paste0("is.DAG(A_succ)")
      val = eval(parse(text = act_to_eval))
    
    if (val != 0){
      break
    }
  }
  
  A_new = A_succ
  
  return(list(A_new = A_new, type.operator = O[i,1], nodes = O[i,2:3]))
  
}

# 3 actions

# A     adjacency matrix of CPDAG C
# x,y   nodes involved in the action


id = function(A, nodes){
  x = nodes[1]
  y = nodes[2]
  A[x,y] = 1
  return(A)
}

dd = function(A, nodes){
  x = nodes[1]
  y = nodes[2]
  A[x,y] = 0
  return(A)
}

rd = function(A, nodes){ # reverse D x -> y
  x = nodes[1]
  y = nodes[2]
  A[x,y] = 0 
  A[y,x] = 1
  return(A)
}
