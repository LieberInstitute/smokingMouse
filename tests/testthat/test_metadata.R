context("metadata validity")

test_that("metadata is valid",
{
    pathToPackage = "/Users/daiannagonzalez/Desktop/smokingMouse"
    expect_type(AnnotationHubData::makeAnnotationHubMetadata(pathToPackage , "metadata.csv"), "list")
})
