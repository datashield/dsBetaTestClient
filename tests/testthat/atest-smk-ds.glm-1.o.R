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

context("dsBetaTestClient::ds.glm.o 1:smoke")

source("connection_to_datasets/init_all_datasets.R")
source("connection_to_datasets/init_smk_datasets.R")

connect.smk.dataset.survival(list("survtime", "time.id", "female", "age.60"))

#
# Tests
#

context("dsBetaTestClient::ds.glm.o(): Standard Poisson regression model for piecewise exponential regression analysis:smoke")

# mod.D<-ds.glm.o("D$survtime~D$time.id+D$female+D$age.60",family="poisson")
mod.D<-ds.glm.o("D$survtime~1+D$time.id+D$female",family="poisson")

print("----")
print(mod.D)
print("----")

output.D<-c(mod.D$coefficients[,1],mod.D$coefficients[,2])

test_that("glm_poisson", {
    expect_equal(ds.ls()$sim1[2],output.D,output.R)
})

context("dsBetaTestClient::ds.glm.o() errors:smoke")
test_that("glm_errors", {
    expect_error(ds.glm.o(), "argument is of length zero", fixed=TRUE)
})

#
# Done
#

context("dsBetaTestClient::ds.glm.o 1:smoke done")
