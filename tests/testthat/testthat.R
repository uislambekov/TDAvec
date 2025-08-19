options(rgl.useNULL=TRUE)

library(testthat)
library(TDAvec)
library(TDAstats)

load_unit_circle <- function() {
  file_path <- system.file("extdata", "unitCircle.csv", package = "TDAvec")
  data <- read.csv(file_path)
  return(data);
}


cloud <- load_unit_circle()
diag = TDAstats::calculate_homology(cloud, threshold = 2)
diag <- rbind(diag, c(0,0,2))

scaleSeq = seq(from=0, to=2, length.out = 11)

test_that("PL_0", {
  R <- as.vector(computePersistenceLandscape(diag, 0, scaleSeq))
  python  <- c( 0, 0.2, 0.4, 0.6, 0.8, 1, 0.8, 0.6, 0.4, 0.2, 0)
  expect_equal(R, python)
})

test_that("PL_1", {
  R <- as.vector(computePersistenceLandscape(diag, 1, scaleSeq))
  python  <- c( 0, 0.0142176303784579, 0.000931983027759931, 0.191114535909928,
                0.391114535909928, 0.269212150161859, 0.0692121501618592, 0, 0, 0, 0)
  expect_equal(R, python)
})


test_that("PS_0", {
R <- computePersistenceSilhouette(diag, 0, scaleSeq)
  python  <- c(0.0507014405880279, 0.0420473568616778, 0.0615450397748265, 0.0861630556847571,
           0.110781071594688, 0.110781071594688, 0.0861630556847571,
           0.0615450397748265, 0.0369270238648959, 0.0123090079549653)
  expect_equal(R, python)
})



test_that("PS_1", {
  R <- computePersistenceSilhouette(diag, 1, scaleSeq)
  python  <- c(
    0.000119825669467705, 0.00125047122181634, 0.0651147889829073, 0.207585992499712,
    0.257838809390676, 0.120660660329335, 0.00853962588653096, 0, 0, 0    )
  expect_equal(R, python)
})

test_that("NL_0", {
  R <- computeNormalizedLife(diag, 0, scaleSeq)
  python  <- c(
    0.817130850366702, 0.234100823784605, 0.123090079549653, 0.123090079549653,
    0.123090079549653, 0.123090079549653, 0.123090079549653, 0.123090079549653,
    0.123090079549653, 0.123090079549653
  )
  expect_equal(R, python)
})

test_that("NL_1", {
  R <- computeNormalizedLife(diag, 1, scaleSeq)
  python  <- c(
    0.0130993953929026, 0.0631855672482697, 0.682417578522476,
    0.713073264620286, 0.713073264620286, 0.713073264620286, 0.246766669336532,
    0, 0, 0  )
  expect_equal(R, python)
})

test_that("VAB_0", {
  R <- computeBettiCurve(diag, 0, scaleSeq)
  python  <- c(
    65.9281367702083, 7.31317883042165, 1, 1, 1, 1, 1, 1, 1, 1  )
  expect_equal(R, python)
})

test_that("VAB_1", {
  R <- computeBettiCurve(diag, 1, scaleSeq)
  python  <- c(
    0.284203374082668, 1.37424398158854, 1.02801848480822, 1, 1, 1, 0.346060750809296, 0, 0, 0  )
  expect_equal(R, python)
})

test_that("ECC", {
  R <- computeEulerCharacteristic(diag, scaleSeq)
  python  <- c(
    65.6439333961256, 5.93893484883312, -0.0280184848082246, 0, 0, 0, 0.653939249190704, 1, 1, 1
    )
  expect_equal(R, python)
})

test_that("PES_0", {
  R <- computePersistentEntropy(diag, 0, scaleSeq)
  python  <- c(
    4.82683006650221, 1.01731577228001, 0.372004512741568, 0.372004512741568, 0.372004512741568,
    0.372004512741568, 0.372004512741568, 0.372004512741568, 0.372004512741568, 0.372004512741568)
  expect_equal(R, python)
})

test_that("PES_1", {
  R <- computePersistentEntropy(diag, 1, scaleSeq)
  python  <- c(
    0.0567767350968682, 0.264582251401072, 0.338458609058045, 0.347892602095769, 0.347892602095769,
    0.347892602095769, 0.120391975082262, 0, 0, 0)
  expect_equal(R, python)
})


test_that("VPB_0", {
  diag_ <- diag
  diag_[,3] <- diag_[,3] - diag_[,2]
  ySeqH0 <- quantile(diag_[diag_[,1]==0, 3], seq(0, 1.1, 0.2))
  R <- computePersistenceBlock(diag_, 0, xSeq=NA, ySeq = ySeqH0)
  python  <- c(
    0.809001048929071, 3.64481540058609, 5.48661206724003, 6.38281631026758, 1.01395250734014)
  expect_equal(R, python)
})


test_that("VPB_1", {
  diag_ <- diag
  diag_[,3] <- diag_[,3] - diag_[,2]
  xSeqH1 <- quantile(diag_[diag_[,1]==1, 2], seq(0, 1.1, 0.2))
  ySeqH1 <- quantile(diag_[diag_[,1]==1, 3], seq(0, 1.1, 0.2))
  R <- computePersistenceBlock(diag_, 1, xSeq=xSeqH1, ySeq = ySeqH1)
  python  <- c(
    0, 0.000499559782281428, 0.00903377330326861, 0.0305098817964819,
    0.00579736873795015, 0.00951944957273334, 0.0261053156675721, 0.0203673036116531,
    0.00428701302896287, 0.0107446149362357, 0.0425768331028477, 0.105346832126612,
    0.145880724120365, 0, 0, 0.0328035481022359, 0.0276785592055745, 0.108908859370022,
    0.131882225448411, 0.0168304894683267, 0.293939456124958, 0.316298612637143, 0.334451098933436,
    0.356748693217673, 0.36362269411629)
  expect_equal(R, python)
})
