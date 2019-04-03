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

context("dsBetaTestClient::ds.asFactor.o")

logindata <- DSLite::setupSURVIVALTest("dsBetaTest", env = environment())
options(datashield.variables=list('survtime', 'time.id', 'female', 'age.60'))
conns <- datashield.login(logins=logindata, assign=TRUE, variables=getOption("datashield.variables", NULL))

#
# Tests
#

ds.asNumeric("D$time.id","TID")

context("dsBetaTestClient::ds.asFactor.o(force.factor.levels)")

test_that("with no force.factor.levels", {
    ds.asFactor.o("TID", "TID.f")

    expect_equal("factor", ds.class("TID.f")$`survival1`)
    expect_equal("factor", ds.class("TID.f")$`survival2`)
    expect_equal("factor", ds.class("TID.f")$`survival3`)
    expect_equal("All tables are valid!", ds.table1D("TID.f")$validity)
})

test_that("with forced.factor.levels of 1:6", {
    ds.asFactor.o("TID", "TID.f2", forced.factor.levels=1:6)

    expect_equal("factor", ds.class("TID.f2")$`survival1`)
    expect_equal("factor", ds.class("TID.f2")$`survival2`)
    expect_equal("factor", ds.class("TID.f2")$`survival3`)
    expect_equal("All tables are valid!", ds.table1D("TID.f2")$validity)
})

test_that("with force.factor.levels of 0:10", {
    ds.asFactor.o("TID", "TID.f3", forced.factor.levels=0:10)

    expect_equal("factor", ds.class("TID.f3")$`survival1`)
    expect_equal("factor", ds.class("TID.f3")$`survival2`)
    expect_equal("factor", ds.class("TID.f3")$`survival3`)
    expect_equal("All tables are valid!", ds.table1D("TID.f3")$validity)
})

test_that("with force.factor.levels of 2:3", {
    ds.asFactor.o("TID", "TID.f4", forced.factor.levels=2:3)

    expect_equal("factor", ds.class("TID.f4")$`survival1`)
    expect_equal("factor", ds.class("TID.f4")$`survival2`)
    expect_equal("factor", ds.class("TID.f4")$`survival3`)
    expect_equal("All tables are valid!", ds.table1D("TID.f4")$validity)
})

test_that("with force.factor.levels of c(1,2,3,4,'a','h',5)", {
    ds.asFactor.o("TID", "TID.f5", forced.factor.levels=c(1,2,3,4,'a','h',5))

    expect_equal("factor", ds.class("TID.f5")$`survival1`)
    expect_equal("factor", ds.class("TID.f5")$`survival2`)
    expect_equal("factor", ds.class("TID.f5")$`survival3`)
    expect_equal("All tables are valid!", ds.table1D("TID.f5")$validity)
})

context("dsBetaTestClient::ds.asFactor.o(fixed.dummy.vars)")

test_that("with fixed.dummy.vars of TRUE", {
    ds.asFactor.o("TID", "TID.mat1", fixed.dummy.vars=TRUE)

    expect_equal("matrix", ds.class("TID.mat1")$`survival1`)
    expect_equal("matrix", ds.class("TID.mat1")$`survival2`)
    expect_equal("matrix", ds.class("TID.mat1")$`survival3`)
})

test_that("with fixed.dummy.vars of TRUE and baseline.level of 6", {
    ds.asFactor.o("TID", "TID.mat6", fixed.dummy.vars=TRUE,baseline.level=6)

    expect_equal("matrix", ds.class("TID.mat6")$`survival1`)
    expect_equal("matrix", ds.class("TID.mat6")$`survival2`)
    expect_equal("matrix", ds.class("TID.mat6")$`survival3`)
})

#
# Tear down
#

datashield.logout(conns)

context("dsBetaTestClient::ds.asFactor.o done")
