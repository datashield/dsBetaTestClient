source("connection_to_datasets/init_all_datasets.R")


context("Returns expected numerical value")
test_that("combined data set",
{
    connect.all.datasets()
    factor.from.file <- list(levels(as.factor(ds.test_env$local.values[,15])))
    factor.from.server <- ds.asFactor.o('D$FACTOR_INTEGER','FACTOR_INT.f', datasources = ds.test_env$connection.opal)
#    print(class(factor.from.server[1]))
#    print(class(factor.from.file))
    expect_true(length(setdiff(factor.from.file, factor.from.server[1]))==0)
})

test_that("split dataset",
{
    connect.dataset.1()
    factor.from.file <- list(levels(as.factor(ds.test_env$local.values.1[,15])))
    factor.from.server <- ds.asFactor.o('D$FACTOR_INTEGER','FACTOR_INT.f', datasources = ds.test_env$connection.opal)
#    print(class(factor.from.server[1]))
#    print(class(factor.from.file))
    expect_true(length(setdiff(factor.from.file, factor.from.server[1]))==0)
    
    connect.dataset.2()
    factor.from.file <- list(levels(as.factor(ds.test_env$local.values.2[,15])))
    factor.from.server <- ds.asFactor.o('D$FACTOR_INTEGER','FACTOR_INT.f', datasources = ds.test_env$connection.opal)
    expect_true(length(setdiff(factor.from.file, factor.from.server[1]))==0)
    
    connect.dataset.3()
    factor.from.file <- list(levels(as.factor(ds.test_env$local.values.3[,15])))
    factor.from.server <- ds.asFactor.o('D$FACTOR_INTEGER','FACTOR_INT.f', datasources = ds.test_env$connection.opal)
#    print(class(factor.from.server[1]))
#    print(class(factor.from.file))
    expect_true(length(setdiff(factor.from.file, factor.from.server[1]))==0)
})


context("Mathematical properties")
test_that("every factor should be unique",
{
  connect.all.datasets()
  factor.from.server <- ds.asFactor.o('D$FACTOR_INTEGER','FACTOR_INT.f', datasources = ds.test_env$connection.opal)
  factors.vector <- unlist(factor.from.server[1])
  expect_true(!any(duplicated(factors.vector)))
})



