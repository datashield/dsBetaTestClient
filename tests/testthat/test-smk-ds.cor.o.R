#-------------------------------------------------------------------------------
# Copyright (c) 2014 OBiBa,
#               2018 University of Newcastle upon Tyne. All rights reserved.
#
# This program and the accompanying materials
# are made available under the terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------------

#
# Set up
#

context("dsBetaTestClient::ds.cor.o:smoke")

source("connection_to_datasets/init_all_datasets.R")
source("connection_to_datasets/init_smk_datasets.R")

connect.smk.dataset.survival(list("survtime", "time.id", "female"))

#
# Tests
#

context("dsBetaTestClient::ds.cor.o():smoke")

test_that("simple test", {
    res <- ds.cor.o(x="D$survtime", y="D$time.id")

    expect_length(res, 3)
    expect_length(res[[1]], 5)
    expect_equal(class(res[[1]]$`Number of missing values in each variable`), "matrix")
    expect_equal(class(res[[1]]$`Number of missing values pairwise`), "matrix")
    expect_equal(class(res[[1]]$`Correlation Matrix`), "matrix")
    expect_equal(class(res[[1]]$`Number of complete cases used`), "matrix")
    expect_equal(res[[1]]$`Error message`, NA)
    expect_length(res[[2]], 5)
    expect_equal(class(res[[2]]$`Number of missing values in each variable`), "matrix")
    expect_equal(class(res[[2]]$`Number of missing values pairwise`), "matrix")
    expect_equal(class(res[[2]]$`Correlation Matrix`), "matrix")
    expect_equal(class(res[[2]]$`Number of complete cases used`), "matrix")
    expect_equal(res[[2]]$`Error message`, NA)
    expect_length(res[[3]], 5)
    expect_equal(class(res[[3]]$`Number of missing values in each variable`), "matrix")
    expect_equal(class(res[[3]]$`Number of missing values pairwise`), "matrix")
    expect_equal(class(res[[3]]$`Correlation Matrix`), "matrix")
    expect_equal(class(res[[3]]$`Number of complete cases used`), "matrix")
    expect_equal(res[[3]]$`Error message`, NA)

})

#
# Done
#

context("dsBetaTestClient::ds.cor.o 1:smoke done")
