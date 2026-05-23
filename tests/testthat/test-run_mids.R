test_that("single observation returns one individual", {
  
  x <- data.frame(
    CameraID="C1",
    Species="A",
    Frame=1,
    Length=100,
    Precision_Error=2
  )
  
  result <- run_mids(
    data=x,
    group_cols=c("CameraID","Species"),
    frame_col="Frame",
    length_col="Length",
    precision_col="Precision_Error"
  )
  
  expect_equal(nrow(result),1)
  
})

test_that("distinct fish remain distinct", {
  
  x <- data.frame(
    CameraID=c("C1","C1"),
    Species=c("A","A"),
    Frame=c(1,2),
    Length=c(100,150),
    Precision_Error=c(2,2)
  )
  
  result <- run_mids(
    data=x,
    group_cols=c("CameraID","Species"),
    frame_col="Frame",
    length_col="Length",
    precision_col="Precision_Error"
  )
  
  expect_equal(nrow(result),2)
  
})

test_that("repeat observations collapse", {
  
  x <- data.frame(
    CameraID=c("C1","C1"),
    Species=c("A","A"),
    Frame=c(1,2),
    Length=c(100,101),
    Precision_Error=c(3,3)
  )
  
  result <- run_mids(
    data=x,
    group_cols=c("CameraID","Species"),
    frame_col="Frame",
    length_col="Length",
    precision_col="Precision_Error"
  )
  
  expect_equal(nrow(result),1)
  
})