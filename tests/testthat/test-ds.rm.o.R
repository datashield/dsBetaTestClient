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

context("dsBetaTestClient::ds.rm.o")

logindata <- DSLite::setupCNSIMTest("dsBetaTest", env = environment())
options(datashield.variables=list("LAB_TSC"))
conns <- datashield.login(logins=logindata, assign=TRUE, variables=getOption("datashield.variables", NULL))

#
# Tests
#

context("dsBetaTestClient::ds.rm.o()")
test_that("mean values [combine]", {
    rm.list <- ds.rm.o('D$LAB_TSC')

    expect_true(length(rm.list) != 0)
})

#
# Tear down
#

datashield.logout(conns)

context("dsBetaTestClient::ds.rm.o done")
