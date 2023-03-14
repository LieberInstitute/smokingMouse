context("metadata validity")

test_that("metadata is valid",
{
    path <- find.package("smokingMouse")
    metadata <- "metadata.csv"
    expect_type(AnnotationHubData::makeAnnotationHubMetadata(path, metadata), "list")
})
