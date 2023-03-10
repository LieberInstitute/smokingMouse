context("metadata validity")

test_that("metadata is valid",
{
    if(!requireNamespace("AnnotationHubData", quietly = TRUE))
        BiocManager::install("AnnotationHubData")

    path <- find.package("smokingMouse")
    metadata <- "metadata.csv"
    expect_type(AnnotationHubData::makeAnnotationHubMetadata(path, metadata), "list")
})
