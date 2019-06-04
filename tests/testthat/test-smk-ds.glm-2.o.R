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

context("dsBetaTestClient::ds.glm.o 2:smoke")

source("connection_to_datasets/init_all_datasets.R")
source("connection_to_datasets/init_smk_datasets.R")

connect.smk.dataset.sim(list("LAB_TSC", "LAB_TRIG"))

#
# Tests
#

context("dsBetaTestClient::ds.glm.o(): Standard Gaussian regression model for piecewise exponential regression analysis:smoke")

mod.D <- ds.glm.o('D$LAB_TSC~D$LAB_TRIG',family="gaussian")
output.D1 <- mod.D$coefficients[,1]
output.D2 <- mod.D$coefficients[,2]

output.R1 <- c('(Intercept)' = 5.697114181, 'D$LAB_TRIG' = 0.075299766)
output.R2 <- c('(Intercept)' = 0.019983645, 'D$LAB_TRIG' = 0.007754159)

test_that("glm_gaussian", {
    expect_equal(output.D1, output.R1)
    expect_equal(output.D2, output.R2)
})

#
# Done
#

context("dsBetaTestClient::ds.glm.o:smoke done")
