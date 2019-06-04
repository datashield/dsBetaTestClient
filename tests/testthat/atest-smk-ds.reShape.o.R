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

context("dsBetaTestClient::ds.reShape.o:smoke")

source("connection_to_datasets/init_all_datasets.R")
source("connection_to_datasets/init_smk_datasets.R")

connect.smk.dataset.sim(list("LAB_TSC", "LAB_HDL"))

#
# Tests
#

context("dsBetaTestClient::ds.reShape.o()")
test_that("simplest ds.reShape.o", {
    myvectors <- c("D$LAB_TSC", "D$LAB_HDL")
    ds.dataFrame.o(x=myvectors, newobj = "orig_df")

    reshape.res <- ds.reShape.o("orig_df", varying=c("A"), v.names=c("B"), direction='wide')

    print(paste0("[" + reshape.res + "]"))

    expect_true(length(reshape.res) == 0)
})

#
# Done
#

context("dsBetaTestClient::ds.reShape.o:smoke done")
