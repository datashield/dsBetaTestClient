#-------------------------------------------------------------------------------
# Copyright (c) 2014 OBiBa,
#               2019 University of Newcastle upon Tyne. All rights reserved.
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

context("dsBetaTestClient::ds.dim.o:smoke")

source("connection_to_datasets/init_all_datasets.R")
source("connection_to_datasets/init_smk_datasets.R")

connect.smk.dataset.sim(list("LAB_TSC"))

#
# Tests
#

context("dsBetaTestClient::ds.dim.o():smoke")
test_that("simple dim", {
    dim.res <- ds.dim.o('D')

    expect_equal(dim.res$`dimensions of D in sim1`[[1]], 2163)
    expect_equal(dim.res$`dimensions of D in sim2`[[1]], 3088)
    expect_equal(dim.res$`dimensions of D in sim3`[[1]], 4128)
    expect_equal(dim.res$`dimensions of D in combined studies`[[1]], 9379)
})

#
# Done
#

context("dsBetaTestClient::ds.dim.o:smoke done")
