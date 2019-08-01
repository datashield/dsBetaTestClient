
.test.basic.expectation <- function(size,variable.created)
{
   if (size <= 0)
   {
     expect_error(ds.rPois.o(samp.size = size,newobj=variable.created, datasources=ds.test_env$connection.opal))
   }
   else if (size > (2^31-1))
   {
     expect_error(ds.rPois.o(samp.size = size,newobj=variable.created, datasources=ds.test_env$connection.opal))
   }
   else
   {
   
      ds.rPois.o(samp.size = size,newobj=variable.created, datasources=ds.test_env$connection.opal)
      exists <- ds.exists(variable.created,datasources=ds.test_env$connection.opal)
      type <- ds.class(variable.created,datasources=ds.test_env$connection.opal)
      sample.size <- ds.length(variable.created, type = "split", datasources=ds.test_env$connection.opal)

     
      for(index in 1:length(exists))
      {
        expect_true(exists[[index]])
      }
   
      for(index in 1:length(sample.size))
      {
        expect_true(sample.size[[index]] == size)
      }
   
      expect_true(type[[1]][1]=="integer")
   }
}

.test.basic.expectation.with.seeds <- function(size,variable.created)
{
  seed <- as.integer(as.POSIXct(Sys.time(), "GMT"))/1000
 
  if (size <= 0)
  {
    expect_error(ds.rPois.o(samp.size = size, 
                            newobj=variable.created,
                            seed.as.integer = seed,
                            return.full.seed.as.set = FALSE,
                            datasources=ds.test_env$connection.opal))
  }
  else if (size > (2^31-1))
  {
    expect_error(ds.rPois.o(samp.size = size,
                            newobj=variable.created, 
                            seed.as.integer = seed,
                            return.full.seed.as.set = FALSE,
                            datasources=ds.test_env$connection.opal))
  }
  
  else
  {
  
    ds.rPois.o(samp.size = size, 
               newobj=variable.created,
               seed.as.integer = seed,
               return.full.seed.as.set = FALSE,
               datasources=ds.test_env$connection.opal)
    exists <- ds.exists(variable.created,datasources=ds.test_env$connection.opal)
    type <- ds.class(variable.created,datasources=ds.test_env$connection.opal)
    sample.size <- ds.length(variable.created, type = "split", datasources=ds.test_env$connection.opal)
    
    for(index in 1:length(exists))
    {
      expect_true(exists[[index]])
    }
    
    for(index in 1:length(sample.size))
    {
      expect_true(sample.size[[index]] == size)
    }
    
    expect_true(type[[1]][1]=="integer")
  }
 
}

.test.too.large.seed <- function(no.server)
{
  seed <- .Machine$integer.max
  if (seed > seed/(no.server +1))
  {
    expect_error(ds.rPois.o(samp.size = size,
                            newobj=variable.created, 
                            seed.as.integer = seed,
                            return.full.seed.as.set = FALSE,
                            datasources=ds.test_env$connection.opal))
  }
}

.test.lambdas <- function(size, variable.created, lambdas)
{
  seed <- as.integer(as.POSIXct(Sys.time(), "GMT"))/1000

  ds.rPois.o(samp.size = size, lambda = c(10,10,10),  seed.as.integer = seed,newobj=variable.created, datasources=ds.test_env$connection.opal)
 
  sample.size <- ds.length(variable.created, type = "split", datasources=ds.test_env$connection.opal)
    
  for(index in 1:length(sample.size))
  {
    expect_true(sample.size[[index]] == lamgdas[i])
  }
}