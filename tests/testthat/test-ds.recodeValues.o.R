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

context("dsBetaTestClient::ds.recodeValues.o")

logindata <- DSLite::setupSURVIVALTest("dsBetaTest", env = environment())
options(datashield.variables=list('survtime', 'time.id', 'female', 'age.60'))
conns <- datashield.login(logins=logindata, assign=TRUE, variables=getOption("datashield.variables", NULL))

#
# Tests
#

context("dsBetaTestClient::ds.recodeValues.o(...) valid arguments")
test_that("recodeValues_checks", {
    res <- ds.recodeValues.o(var.name="D$age.60", values2replace.vector=c(-1,-6), new.values.vector=c(-10,-16), datasources=conns)

    expect_equal(res$is.object.created, "A data object <D$age.60_recoded> has been created in all specified data sources")
    expect_equal(res$validity.check, "<D$age.60_recoded> appears valid in all sources")
})

context("dsBetaTestClient::ds.recodeValues.o(...) invalid arguments")
test_that("recodeValues_erros", {
    expect_error(ds.recodeValues.o(var.name=NULL), "Please provide the name of the variable to be recoded: eg 'xxx'")
    expect_error(ds.recodeValues.o(var.name="age.60", values2replace.vector=NULL), "Please provide a vector specifying the values to be replaced eg c\\(1\\,7\\,NA\\)")
    expect_error(ds.recodeValues.o(var.name="age.60", values2replace.vector=c(-1, -6), new.values.vector=NULL), "Please provide a vector specifying the new values to be set eg c\\(3\\,NA\\,4\\)")
    expect_error(ds.recodeValues.o(var.name="age.60", values2replace.vector=c(-1, -6, -8), new.values.vector=c(-10, -16)), "Please ensure that values2replace.vector and new.values.vector have same length and are in the same order")
    expect_error(ds.recodeValues.o(var.name="age.60", values2replace.vector=c(-1, -6, -6), new.values.vector=c(-10, -16, -17)), "No value may appear more than once in the values2replace.vector")
})

#
# Tear down
#

datashield.logout(conns)

context("dsBetaTestClient::ds.recodeValues.o done")
