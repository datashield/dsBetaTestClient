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

context("dsBetaTestClient::ds.tapply.o:smoke")

source("connection_to_datasets/init_all_datasets.R")
source("connection_to_datasets/init_smk_datasets.R")

connect.smk.dataset.sim(list("LAB_TSC", "GENDER"))

#
# Tests
#

context("dsBetaTestClient::ds.tapply.o()")
test_that("simplest 'ds.tapply.o'", {
    list <- ds.tapply.o('D$LAB_TSC', INDEX.names=c('D$GENDER'), FUN.name='sum', datasources=ds.test_env$connection.opal)

    expect_true(length(list) == 0)
})

#
# Tear down
#

context("dsBetaTestClient::ds.tapply.o:smoke done")
