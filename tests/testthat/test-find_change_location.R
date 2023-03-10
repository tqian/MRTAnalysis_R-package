# generate a vector containing different values
v <- c("a", "a", "b", "c", "c")

test_that("multiplication works", {
    expect_equal(find_change_location(v), c(1, 3, 4))
})

v_test <- {}
test_that(
    "check error for find_change_location",
    {
        expect_error(
            find_change_location(v_test),
            "The vector need to have length > 1."
        )
    }
)
