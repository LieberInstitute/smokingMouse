context("metadata validity")

test_that("metadata is valid",
{
    if(!requireNamespace("AnnotationHubData", quietly = TRUE))
        BiocManager::install("AnnotationHubData")

    path <- find.package("smokingMouse")
    metadata <- system.file("extdata", "metadata.csv", package = "smokingMouse")
    expect_true(AnnotationHubData::makeAnnotationHubMetadata(path, metadata))
})
