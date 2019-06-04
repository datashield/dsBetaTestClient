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

context("dsBetaTestClient::ds.mean.o:smoke")

source("connection_to_datasets/init_all_datasets.R")
source("connection_to_datasets/init_smk_datasets.R")

connect.smk.dataset.sim(list("LAB_TSC"))

#
# Tests
#

context("dsBetaTestClient::ds.mean.o:(type=combine):smoke")
test_that("mean values [combine]", {
    stat.mean <- ds.mean.o(x='D$LAB_TSC',type='combine')

    expect_false(is.na(stat.mean$Global.Mean[1]))
    expect_equal(as.numeric(stat.mean$Global.Mean[1]), 5.85192485623003, tolerance = .000000000000001)
})

context("dsBetaTestClient::ds.mean.o(type=combine) loose:smoke")
test_that("mean values [combine] loose", {
    stat.mean <- ds.mean.o(x='D$LAB_TSC',type='combine')

    expect_false(is.na(stat.mean$Global.Mean[1]))
    expect_equal(as.numeric(stat.mean$Global.Mean[1]), 5.85192485623003, tolerance = .000000000000001)
})

context("dsBetaTestClient::ds.mean.o(type=split):smoke")
test_that("mean values [split]", {
    stat.mean <- ds.mean.o(datasources=ds.test_env$connection.opal, x='D$LAB_TSC', type='split')

    expect_false(is.na(stat.mean$Mean.by.Study[1]))
    expect_equal(stat.mean$Mean.by.Study[1], 5.87211344770338, tolerance = .000000000000001)
    expect_false(is.na(stat.mean$Mean.by.Study[2]))
    expect_equal(stat.mean$Mean.by.Study[2], 5.84526388341867, tolerance = .000000000000001)
    expect_false(is.na(stat.mean$Mean.by.Study[3]))
    expect_equal(stat.mean$Mean.by.Study[3], 5.84630008623168, tolerance = .000000000000001)
})

context("dsBetaTestClient::ds.mean.o() test errors:smoke")
test_that("mean_erros", {
    ds.asCharacter.o(x='D$LAB_TSC', newobj="not_a_numeric")

    expect_error(ds.mean.o(), "Please provide the name of the input vector!", fixed=TRUE)
#    expect_error(ds.mean.o(x='D$LAB_TSC', type='datashield'), 'Function argument "type" has to be either "combine" or "split"', fixed=TRUE)
#    expect_error(ds.mean.o(x='not_a_numeric'), "The input object must be an integer or a numeric vector.", fixed=TRUE)
})

#
# Done
#

context("dsBetaTestClient::ds.mean.o:smoke done")
