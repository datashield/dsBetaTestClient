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

context("dsBetaTestClient::ds.glm.o 1")

logindata <- DSLite::setupSURVIVALTest("dsBetaTest", env = environment())
options(datashield.variables=list('survtime', 'time.id', 'female', 'age.60'))
conns <- datashield.login(logins=logindata, assign=TRUE, variables=getOption("datashield.variables", NULL))

#
# Tests
#

context("dsBetaTestClient::ds.glm.o(): Standard Poisson regression model for piecewise exponential regression analysis")
test_that("glm_poisson", {
     mod.D<-ds.glm.o("D$survtime~1+D$time.id+D$female+D$age.60", family="poisson")

    expect_equal(mod.D$Nmissing[1], 134)
    expect_equal(mod.D$Ntotal, 6388)
#    expect_equal(mod.D$errorMessage$sim1, "No errors")
#    expect_equal(mod.D$errorMessage$sim2, "No errors")
#    expect_equal(mod.D$errorMessage$sim3, "No errors")
    expect_equal(mod.D$dev, 3439.32, tolerance = .0005)
})

context("dsBetaTestClient::ds.glm.o() errors")
test_that("glm_errors", {
    expect_error(ds.glm.o(), "Please provide a valid regression formula!", fixed=TRUE)
})

#
# Tear down
#

datashield.logout(conns)

context("dsBetaTestClient::ds.glm.o 1 done")
