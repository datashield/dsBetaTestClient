source("connection_to_datasets/init_all_datasets.R")
connect.all.datasets()

context("Returns expected numerical value")
test_that("combined data set",
{
  #  print("positive and negative integers value")
    mean.from.files <- mean(ds.test_env$local.values[,6])
    mean.from.servers <- ds.mean.o(x='D$INTEGER',type='combine', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)
    expect_equal(mean.from.servers[[1]][1], mean.from.files, tolerance = 10^-15)
   
   # print("positive integers")
    mean.from.files <- mean(ds.test_env$local.values[,8])
    mean.from.servers <- ds.mean.o(x='D$POSITIVE_INTEGER',type='combine', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)
    expect_equal(mean.from.servers[[1]][1], mean.from.files, tolerance = 10^-15)
    
  #  print("negative integers")
    mean.from.files <- mean(ds.test_env$local.values[,9])
    mean.from.servers <- ds.mean.o(x='D$NEGATIVE_INTEGER',type='combine', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)
    expect_equal(mean.from.servers[[1]][1], mean.from.files, tolerance = 10^-15)
    
   # print("non negative integers")
    mean.from.files <- mean(ds.test_env$local.values[,7])
    mean.from.servers <- ds.mean.o(x='D$NON_NEGATIVE_INTEGER',type='combine', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)
    expect_equal(mean.from.servers[[1]][1], mean.from.files, tolerance = 10^-15)
    #----------------------
    #print("positive and negative numeric value")
    mean.from.files <- mean(ds.test_env$local.values[,10])
    mean.from.servers <- ds.mean.o(x='D$NUMERIC',type='combine', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)
    expect_equal(mean.from.servers[[1]][1], mean.from.files, tolerance = 10^-15)
    
    #print("positive numeric")
    mean.from.files <- mean(ds.test_env$local.values[,12])
    mean.from.servers <- ds.mean.o(x='D$POSITIVE_NUMERIC',type='combine', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)
    expect_equal(mean.from.servers[[1]][1], mean.from.files, tolerance = 10^-15)
    
  #  print("negative numeric")
    mean.from.files <- mean(ds.test_env$local.values[,13])
    mean.from.servers <- ds.mean.o(x='D$NEGATIVE_NUMERIC',type='combine', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)
    expect_equal(mean.from.servers[[1]][1], mean.from.files, tolerance = 10^-15)
    
   # print("non negative numeric")
    mean.from.files <- mean(ds.test_env$local.values[,11])
    mean.from.servers <- ds.mean.o(x='D$NON_NEGATIVE_NUMERIC',type='combine', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)
    expect_equal(mean.from.servers[[1]][1], mean.from.files, tolerance = 10^-15)
})

test_that("split dataset",
{
    #print("positive and negative integers value")
    mean.from.file.1 <- mean(ds.test_env$local.values.1[,6])
    mean.from.file.2 <- mean(ds.test_env$local.values.2[,6])
    mean.from.file.3 <- mean(ds.test_env$local.values.3[,6])
    mean.from.servers <- ds.mean.o(x='D$INTEGER',type='split', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)      
    expect_equal(mean.from.servers[[1]][1], mean.from.file.1, tolerance = 10^-15)
    expect_equal(mean.from.servers[[1]][2], mean.from.file.2, tolerance = 10^-15)
    expect_equal(mean.from.servers[[1]][3], mean.from.file.3, tolerance = 10^-15)
    
    #print("positive integers value")
    mean.from.file.1 <- mean(ds.test_env$local.values.1[,8])
    mean.from.file.2 <- mean(ds.test_env$local.values.2[,8])
    mean.from.file.3 <- mean(ds.test_env$local.values.3[,8])
    mean.from.servers <- ds.mean.o(x='D$POSITIVE_INTEGER',type='split', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)      
    expect_equal(mean.from.servers[[1]][1], mean.from.file.1, tolerance = 10^-15)
    expect_equal(mean.from.servers[[1]][2], mean.from.file.2, tolerance = 10^-15)
    expect_equal(mean.from.servers[[1]][3], mean.from.file.3, tolerance = 10^-15)
    
    #print("negative integers")
    mean.from.file.1 <- mean(ds.test_env$local.values.1[,9])
    mean.from.file.2 <- mean(ds.test_env$local.values.2[,9])
    mean.from.file.3 <- mean(ds.test_env$local.values.3[,9])
    mean.from.servers <- ds.mean.o(x='D$NEGATIVE_INTEGER',type='split', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)      
    expect_equal(mean.from.servers[[1]][1], mean.from.file.1, tolerance = 10^-15)
    expect_equal(mean.from.servers[[1]][2], mean.from.file.2, tolerance = 10^-15)
    expect_equal(mean.from.servers[[1]][3], mean.from.file.3, tolerance = 10^-15)
    
    #print("non negative integers")
    mean.from.file.1 <- mean(ds.test_env$local.values.1[,7])
    mean.from.file.2 <- mean(ds.test_env$local.values.2[,7])
    mean.from.file.3 <- mean(ds.test_env$local.values.3[,7])
    mean.from.servers <- ds.mean.o(x='D$NON_NEGATIVE_INTEGER',type='split', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)      
    
    expect_equal(mean.from.servers[[1]][1], mean.from.file.1, tolerance = 10^-15)
    expect_equal(mean.from.servers[[1]][2], mean.from.file.2, tolerance = 10^-15)
    expect_equal(mean.from.servers[[1]][3], mean.from.file.3, tolerance = 10^-15)
    
    #--------------
    #print("positive and negative numeric value")
    mean.from.file.1 <- mean(ds.test_env$local.values.1[,10])
    mean.from.file.2 <- mean(ds.test_env$local.values.2[,10])
    mean.from.file.3 <- mean(ds.test_env$local.values.3[,10])
    mean.from.servers <- ds.mean.o(x='D$NUMERIC',type='split', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)      
    expect_equal(mean.from.servers[[1]][1], mean.from.file.1, tolerance = 10^-15)
    expect_equal(mean.from.servers[[1]][2], mean.from.file.2, tolerance = 10^-15)
    expect_equal(mean.from.servers[[1]][3], mean.from.file.3, tolerance = 10^-15)
    
    #print("positive numeric value")
    mean.from.file.1 <- mean(ds.test_env$local.values.1[,12])
    mean.from.file.2 <- mean(ds.test_env$local.values.2[,12])
    mean.from.file.3 <- mean(ds.test_env$local.values.3[,12])
    mean.from.servers <- ds.mean.o(x='D$POSITIVE_NUMERIC',type='split', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)      
    expect_equal(mean.from.servers[[1]][1], mean.from.file.1, tolerance = 10^-15)
    expect_equal(mean.from.servers[[1]][2], mean.from.file.2, tolerance = 10^-15)
    expect_equal(mean.from.servers[[1]][3], mean.from.file.3, tolerance = 10^-15)
    
    #print("negative numeric")
    mean.from.file.1 <- mean(ds.test_env$local.values.1[,13])
    mean.from.file.2 <- mean(ds.test_env$local.values.2[,13])
    mean.from.file.3 <- mean(ds.test_env$local.values.3[,13])
    mean.from.servers <- ds.mean.o(x='D$NEGATIVE_NUMERIC',type='split', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)      
    expect_equal(mean.from.servers[[1]][1], mean.from.file.1, tolerance = 10^-15)
    expect_equal(mean.from.servers[[1]][2], mean.from.file.2, tolerance = 10^-15)
    expect_equal(mean.from.servers[[1]][3], mean.from.file.3, tolerance = 10^-15)
    
    #print("non negative numeric")
    mean.from.file.1 <- mean(ds.test_env$local.values.1[,11])
    mean.from.file.2 <- mean(ds.test_env$local.values.2[,11])
    mean.from.file.3 <- mean(ds.test_env$local.values.3[,11])
    mean.from.servers <- ds.mean.o(x='D$NON_NEGATIVE_NUMERIC',type='split', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)      
    
    expect_equal(mean.from.servers[[1]][1], mean.from.file.1, tolerance = 10^-15)
    expect_equal(mean.from.servers[[1]][2], mean.from.file.2, tolerance = 10^-15)
    expect_equal(mean.from.servers[[1]][3], mean.from.file.3, tolerance = 10^-15)
})

#context("Returns expected non-numerical values")
#test_that("character",
#{
#  #print("character")
#  mean.from.servers <- ds.mean.o(x='D$CHARACTER',type='combine', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)
#  expect_equal(is.na(mean.from.servers),TRUE)
#})

#test_that("boolean",
#{          
#  #print("boolean")
#  mean.from.files <- mean(ds.test_env$local.values[,3])
#  mean.from.servers <- ds.mean.o(x='D$BOOLEAN',type='combine', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)
#  print("mean.from.files")
#  expect_equal(mean.from.servers,mean.from.files,tolerance = 10^-15)
#})

context("Mathematical properties")
test_that("residual deviation tends to 0",
{
    mean.from.servers <- ds.mean.o(x='D$INTEGER',type='combine', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)
    expect_equal(sum(ds.test_env$local.values[,6] - mean.from.servers[[1]][1]), 0, tolerance = 10^-8)
    
    mean.from.servers <- ds.mean.o(x='D$POSITIVE_INTEGER',type='combine', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)
    expect_equal(sum(ds.test_env$local.values[,8] - mean.from.servers[[1]][1]), 0, tolerance = 10^-8)
   
    mean.from.servers <- ds.mean.o(x='D$NEGATIVE_INTEGER',type='combine', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)
    expect_equal(sum(ds.test_env$local.values[,9] - mean.from.servers[[1]][1]), 0, tolerance = 10^-8)
    
    mean.from.servers <- ds.mean.o(x='D$NON_NEGATIVE_INTEGER',type='combine', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)
    expect_equal(sum(ds.test_env$local.values[,7] - mean.from.servers[[1]][1]), 0, tolerance = 10^-8)
    #----
    mean.from.servers <- ds.mean.o(x='D$NUMERIC',type='combine', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)
    expect_equal(sum(ds.test_env$local.values[,10] - mean.from.servers[[1]][1]), 0, tolerance = 10^-8)
    
    mean.from.servers <- ds.mean.o(x='D$POSITIVE_NUMERIC',type='combine', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)    
    expect_equal(sum(ds.test_env$local.values[,12] - mean.from.servers[[1]][1]), 0, tolerance = 10^-8)
   
    mean.from.servers <- ds.mean.o(x='D$NEGATIVE_NUMERIC',type='combine', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)
    expect_equal(sum(ds.test_env$local.values[,13] - mean.from.servers[[1]][1]), 0, tolerance = 10^-8)
    
    mean.from.servers <- ds.mean.o(x='D$NON_NEGATIVE_NUMERIC',type='combine', check=TRUE,save.mean.Nvalid=FALSE, datasources=ds.test_env$connection.opal)
    expect_equal(sum(ds.test_env$local.values[,11] - mean.from.servers[[1]][1]), 0, tolerance = 10^-8)
})



#context("incorrect arguments")
#test_that("A suitable value is returned with some incorrect arugements",
#{
#    print(class(ds.test_env$connection.opal))
#    print(ds.mean.o())  
#})


