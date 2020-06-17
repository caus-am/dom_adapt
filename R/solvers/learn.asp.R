learn.asp.available_solvers <- c("clingo")

learn.asp <- function (testConfig, solverConfig,
                       graph_filename_template,
                       currentDir = "./../tmp/",
                       indPath = NULL,
                       background_knowledge_file = NULL,
                       tested_independences,
                       verbose=1,
                       evalDirectCauses = FALSE,
                       parallelise=FALSE) {
  # Function that prepares inputs and calls the ASP solver.
  #
  # Results options (where to store the results and tmp files):
  # currentDir           - the tmp folder where the results are stored.
  # aspSetsFilename      - the serialized version of sets, default: "pipeline.pre.asp",
  # indFilename          - the serialized version of independences, default: "pipeline.ind",
  # outputFilename       - the clingo/foxPSL output file, default: "pipeline.ind.clingo"
  
  if (!solverConfig$solver %in% learn.asp.available_solvers) 
    stop("Solver ", solver, " not available.")
  
  # Encoding sets and write them in a different file (so we can reuse the independences file).
  aspSetsFullPath <- paste(currentDir, graph_filename_template, '.pre.asp', sep='')
  n <- testConfig$n
  
  tic() 
  if (solverConfig$solver == "clingo") {
    writeAspSets(n=n, aspSetsFullPath = aspSetsFullPath, tested_independences=tested_independences)
    build_tree(n=n, aspSetsFullPath = aspSetsFullPath, tested_independences=tested_independences)
  }
  encoding_time<-toc()
  
  # Create the fullpaths for each of the tmp files using the current directory.
  outputFullPath <- paste(currentDir, graph_filename_template, '.ind.clingo', sep='')

  ##############################################################################
  tic()
  if (verbose) cat(" - Solving with config", paste(solverConfig, collapse=","), "\n")

  if (solverConfig$solver == "clingo") {
    clingoCmd <- paste("./../ASP/clingo ", solverConfig$solver_conf, " --const n=", n, sep="")
    clingoInputFiles <- paste(aspSetsFullPath, " ", indPath, " ./../ASP/", solverConfig$encode, " ", background_knowledge_file, sep="")
    
    if (solverConfig$multipleMinima == FALSE){
      cmd <- paste(clingoCmd, clingoInputFiles, "| tee " , outputFullPath)
      cat(cmd)
      system(cmd, ignore.stdout = TRUE)
    } else if (solverConfig$multipleMinima == "short") {
      L <- learn.asp.short(clingoCmd, clingoInputFiles, n=n)              
    } else if (solverConfig$multipleMinima == "iterative" || solverConfig$multipleMinima == "iterativeParallel") {
      L <- learn.asp.iterative(clingoCmd=clingoCmd, clingoInputFiles=clingoInputFiles, n=n, graph_filename_template=graph_filename_template,
                               tested_independences=tested_independences, evalDirectCauses=evalDirectCauses, currentDir=currentDir, parallelise=(solverConfig$multipleMinima == "iterativeParallel"))   
    }
  }
  
  sol_file<- outputFullPath
  solving_time<-toc() 
  
  tic()
  if (solverConfig$solver != "clingo" || solverConfig$multipleMinima == FALSE) {
      L <- parse_asp_solution(n=n, solver=solverConfig$solver, encode=solverConfig$encode, sol_file,  solverConfig$multipleMinima)
  }
  
  parsing_time<-toc()
    
  if (all(is.na(L)) ) {
    L<-list()
    solving_time<-Inf
  }
  
  L$encoding_time<-encoding_time
  L$solving_time<-solving_time
  L$parsing_time<-parsing_time
  
  L
}
