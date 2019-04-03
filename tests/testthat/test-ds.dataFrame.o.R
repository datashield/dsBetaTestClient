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

context("dsBetaTestClient::ds.dataFrame.o")

logindata <- DSLite::setupCNSIMTest("dsBetaTest", env = environment())
options(datashield.variables=NULL)
conns <- datashield.login(logins=logindata, assign=TRUE, variables=getOption("datashield.variables", NULL))

#
# Tests
#

context("dsBetaTestClient::ds.dataFrame.o() create a dataframe")
test_that("dataframe_exists", {
    myvectors <- c('D$LAB_TSC', 'D$LAB_HDL')
    ds.dataFrame.o(x=myvectors)
    res <- ds.ls()

    expect_equal("D", res$sim1[1])
    expect_equal("df_new", res$sim1[2])
    expect_equal("D", res$sim2[1])
    expect_equal("df_new", res$sim2[2])
    expect_equal("D", res$sim3[1])
    expect_equal("df_new", res$sim3[2])
})


context("dsBetaTestClient::ds.dataFrame.o() errors")
test_that("dataframe_errors", {
    expect_error(ds.dataframe(), "Please provide the name of the list that holds the input vectors!", fixed=TRUE)
})

#
# Tear down
#

datashield.logout(conns)

context("dsBetaTestClient::ds.dataFrame.o done")
