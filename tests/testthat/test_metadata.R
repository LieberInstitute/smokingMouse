context("metadata validity")

test_that("metadata is valid",
{
    expect_type(AnnotationHubData::makeAnnotationHubMetadata("metadata.csv"), "list")
})
