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

context("dsBetaTestClient::ds.listServersideFunctions.o:smoke")

source("connection_to_datasets/init_all_datasets.R")
source("connection_to_datasets/init_smk_datasets.R")

connect.smk.dataset.sim(list("LAB_TSC", "LAB_HDL"))

#
# Tests
#

context("dsBetaTestClient::ds.listServersideFunctions.o():smoke check results")
test_that("check results", {
    assign.functions <- factor(c(
      "BooleDS.o", "as.character", "as.null", "as.numeric", "asCharacterDS.o", "asDataMatrixDS.o",
      "asFactorDS", "asFactorDS2.o", "asIntegerDS.o", "asListDS", "asLogicalDS.o", "asMatrixDS",
      "asMatrixDS.o", "asNumericDS.o", "attach", "c", "cDS", "cbind",
      "cbindDS.o", "changeRefGroupDS", "complete.cases", "dataFrameDS.o", "dataFrameSortDS.o", "dataFrameSubsetDS2.o",
      "dataframeDS", "exp", "lexisDS", "lexisDS2.o", "lexisDS3.o", "list",
      "listDS", "log", "matrixDS", "matrixDetDS2", "matrixDiagDS", "matrixDimnamesDS", "matrixInvertDS", 
      "matrixMultDS", "matrixTransposeDS", "mergeDS.o", "rBinomDS.o", "rNormDS.o", "rPoisDS.o",
      "rUnifDS.o", "rbindDS.o", "reShapeDS.o", "recodeLevelsDS", "recodeValuesDS2.o", "rep",
      "replaceNaDS", "rowColCalcDS", "seedDS.o", "seqDS.o", "subsetByClassDS", "subsetDS", "sum",
      "tapplyDS.assign.o", "unlist"
    ))
    aggregate.functions <- factor(c(
      "NROW", "alphaPhiDS", "asFactorDS1.o", "asListDS.o",
      "checkNegValueDS", "class", "colnames", "cor.test",
      "corDS", "covDS", "covDS.o", "dataFrameSubsetDS1.o",
      "densityGridDS", "dim", "dimDS", "dimDS.o",
      "exists", "glmDS1", "glmDS1.o", "glmDS2",
      "glmDS2.o", "glmSLMADS2.o", "histogramDS", "is.character",
      "is.factor", "is.list", "is.null", "is.numeric",
      "isNaDS", "isValidDS", "length", "lengthDS.o",
      "levels", "lexisDS1.o", "listDisclosureSettingsDS.o", "ls",
      "matrixDetDS1", "meanDS", "meanDS.o", "meanSdGpDS", "messageDS.o",
      "namesDS", "numNaDS", "quantileMeanDS", "rangeDS",
      "recodeValuesDS1.o", "rilmDS.b", "rmDS.o", "scatterPlotDS.o",
      "scoreVectDS", "setSeedDS.o", "t.test", "tTestFDS2",
      "table1dDS", "table2dDS", "tapplyDS.o", "testObjExistsDS.o",
      "unListDS.o", "varDS", "varDS.o"
    ))

    res <- ds.listServersideFunctions.o()

    expect_length(res, 2)
    expect_length(res$serverside.assign.functions, 3)
    expect_length(res$serverside.aggregate.functions, 3)
    
    sim1.assign.res    <- res$serverside.assign.functions$sim1
    sim1.aggregate.res <- res$serverside.aggregate.functions$sim1
    sim2.assign.res    <- res$serverside.assign.functions$sim2
    sim2.aggregate.res <- res$serverside.aggregate.functions$sim2
    sim3.assign.res    <- res$serverside.assign.functions$sim3
    sim3.aggregate.res <- res$serverside.aggregate.functions$sim3

    expect_length(sim1.assign.res, 58)
    expect_equal(sim1.assign.res, assign.functions)
    expect_length(sim1.aggregate.res, 60)
    expect_equal(sim1.aggregate.res, aggregate.functions)
    expect_length(sim2.assign.res, 58)
    expect_equal(sim2.assign.res, assign.functions)
    expect_length(sim2.aggregate.res, 60)
    expect_equal(sim2.aggregate.res, aggregate.functions)
    expect_length(sim3.assign.res, 58)
    expect_equal(sim3.assign.res, assign.functions)
    expect_length(sim3.aggregate.res, 60)
    expect_equal(sim3.aggregate.res, aggregate.functions)
})

#
# Done
#

context("dsBetaTestClient::ds.listServersideFunctions.o:smoke done")
