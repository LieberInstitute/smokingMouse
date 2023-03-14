context("metadata validity")

test_that("metadata is valid",
{
    expect_type(AnnotationHubData::makeAnnotationHubMetadata(pathToPackage = "/Users/daiannagonzalez/Desktop/smokingMouse" , "metadata.csv"), "list")
})
