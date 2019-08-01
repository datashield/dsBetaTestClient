source("connection_to_datasets/init_all_datasets.R")
source("definition_tests/def-ds.rPois.R")

context("ds.rPois.o()::expt::no seeds::single")
test_that("Poisson without seeds",
{
  connect.dataset.1()
  .test.basic.expectation(-1,'poisson_dist_1')
  .test.basic.expectation(0,'poisson_dist_2')
  .test.basic.expectation(1,'poisson_dist_3')
  .test.basic.expectation(20,'poisson_dist_4')
  .test.basic.expectation(30,'poisson_dist_5')
  .test.basic.expectation(2^31-1,'poisson_dist_6')
})


context("ds.rPoiss.o()::expt::no seeds::multiple")
test_that("Poisson without seeds",
{
  connect.all.datasets()
  .test.basic.expectation(-1,'poisson_dist_1')
  .test.basic.expectation(0,'poisson_dist_2')
  .test.basic.expectation(1,'poisson_dist_3')
  .test.basic.expectation(20,'poisson_dist_4')
  .test.basic.expectation(30,'poisson_dist_5')
  .test.basic.expectation(2^31-1,'poisson_dist_6')
})



context("ds.rPois.o()::expt::with seeds::single")
test_that("Poisson with seeds",
{
  
  connect.dataset.1()
  
  .test.basic.expectation.with.seeds(-1,'poisson_dist_1')
  .test.basic.expectation.with.seeds(0,'poisson_dist_2')
  .test.basic.expectation.with.seeds(1,'poisson_dist_3')
  .test.basic.expectation.with.seeds(20,'poisson_dist_4')
  .test.basic.expectation.with.seeds(30,'poisson_dist_5')
  .test.basic.expectation.with.seeds(2^31-1,'poisson_dist_6')
  .test.basic.expectation.with.seeds(2^31,'poisson_dist_7')
})

context("ds.rPoiss.o()::expt::with seeds::multiple")
test_that("Poisson with seeds",
{
  connect.all.datasets()
  .test.basic.expectation.with.seeds(-1,'poisson_dist_1')
  .test.basic.expectation.with.seeds(0,'poisson_dist_2')
  .test.basic.expectation.with.seeds(1,'poisson_dist_3')
  .test.basic.expectation.with.seeds(20,'poisson_dist_4')
  .test.basic.expectation.with.seeds(30,'poisson_dist_5')
  .test.basic.expectation.with.seeds(2^31-1,'poisson_dist_6')
  .test.basic.expectation.with.seeds(2^31,'poisson_dist_7')
  .test.too.large.seed(3)
})

context('ds.rPois.o()::ext::lambdas::single')
test_that("lambdas",
{
  connect.all.datasets()
  .test.lambdas(30, 'poisson_dist_1', c(10,10,10))
  .test.lambdas(30, 'poisson_dist_2', c(5,15,10))
})


