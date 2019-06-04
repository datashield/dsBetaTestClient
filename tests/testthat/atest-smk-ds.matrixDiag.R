#-------------------------------------------------------------------------------
# Copyright (c) 2014 OBiBa,
#               2019 University of Newcastle upon Tyne. All rights reserved.
#
# This program and the accompanying materials
# are made available under the terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnurg/licenses/>.
#-------------------------------------------------------------------------------

#
# Set up
#

context("dsBetaTestClient::ds.matrixDiag:smoke")

source("connection_to_datasets/init_all_datasets.R")
source("connection_to_datasets/init_smk_datasets.R")

connect.smk.dataset.sim(list("LAB_TSC"))

#
# Tests
#

context("dsBetaTestClient::ds.matrixDiag()")
test_that("simplest ds.matrixDiag", {
    matrix <- c(0, 1, 2, 3, 4, 5)

    res <- ds.matrixDiag(x1=matrix, aim="clientside.vector.2.matrix", nrows.scalar=6)

    expect_length(res, 2)
    expect_equal(res[[1]], "A data object <new_matrix> has been created in all specified data sources")
    expect_equal(res[[2]], "<new_matrix> appears valid in all sources")

    check.class<-ds.class("new_matrix",datasources=ds.test_env$connection.opal)

    expect_length(check.class, 3)
    expect_equal(check.class$sim1, "matrix")
    expect_equal(check.class$sim2, "matrix")
    expect_equal(check.class$sim3, "matrix")
})

#
# Tear down
#

context("dsBetaTestClient::ds.matrixDiag:smoke done")
