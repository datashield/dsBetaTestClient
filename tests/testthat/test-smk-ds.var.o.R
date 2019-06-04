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

context("dsBetaTestClient::ds.var.o:smoke")

source("connection_to_datasets/init_all_datasets.R")
source("connection_to_datasets/init_smk_datasets.R")

connect.smk.dataset.sim(list("LAB_TSC"))

#
# Tests
#

context("dsBetaTestClient::ds.var.o():smoke")
test_that("simple var", {
    var.res <- ds.var.o(x = 'D$LAB_TSC', datasources = ds.test_env$connection.opal)

    expect_equal(var.res$Variance.by.Study[1], 1.229163, tolerance = .000001)
    expect_equal(var.res$Variance.by.Study[2], 1.140606, tolerance = .000001)
    expect_equal(var.res$Variance.by.Study[3], 1.134995, tolerance = .000001)
})

#
# Done
#

context("dsBetaTestClient::ds.var.o:smoke done")
