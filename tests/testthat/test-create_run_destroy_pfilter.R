
test_that("'create_run_destroy_pfilter' works - no region", {
    set.seed(0)
    for (has_stock_init in c(TRUE, FALSE)) {
        for (is_dominant in c(TRUE, FALSE)) {
            for (cohort_type in c("existing", "new", "expiring")) {
                n_particle <- 5L
                work_dir <- tempfile()
                dir.create(work_dir)
                df <- simulated_df_noreg(
                    has_stock_init = has_stock_init,
                    is_dominant = is_dominant,
                    n_particle = n_particle,
                    is_forecast = FALSE,
                    cohort_type = cohort_type
                )
                create_run_destroy_pfilter(
                    df = df,
                    threshold = 0.5,
                    work_dir = work_dir,
                    is_forecast = FALSE,
                    quiet = TRUE
                )
                filename <- file.path(work_dir, "tmp-population-2000-Female.bin")
                con <- file(filename, "rb")
                n_meta <- readBin(con, what = "integer", n = 1)
                n_count <- readBin(con, what = "integer", n = 1)
                meta <- unserialize(readBin(con, what = "raw", n = n_meta))
                count <- readBin(con, what = "double", n = n_count)
                expect_identical(meta$cohort[[1]], 2000L)
                expect_identical(meta$sexgender[[1]], "Female")
                expect_false(anyNA(count))
                close(con)
                for (nm in c("deaths", "immigration1", "emigration1",
                             "immigration2", "emigration2")) {
                    filename <- file.path(work_dir,
                                          sprintf("tmp-%s-2000-Female.bin", nm))
                    con <- file(filename, "rb")
                    n_meta <- readBin(con, what = "integer", n = 1)
                    n_count <- readBin(con, what = "integer", n = 1)
                    meta <- unserialize(readBin(con, what = "raw", n = n_meta))
                    count <- readBin(con, what = "double", n = n_count)
                    expect_identical(meta$cohort[[1]], 2000L)
                    expect_identical(meta$sexgender[[1]], "Female")
                    expect_false(anyNA(count))
                    close(con)
                }
                diag <- readRDS(file.path(work_dir, "tmp-diagnostics-2000-Female.rds"))
                expect_false(anyNA(diag))
                unlink(work_dir, recursive = TRUE)
            }
        }
    }
})

test_that("'do_forecast' works - no region", {
    set.seed(0)
    for (is_dominant in c(TRUE, FALSE)) {
        for (cohort_type in c("existing", "new", "expiring")) {
            n_particle <- 5L
            work_dir <- tempfile()
            dir.create(work_dir)
            df <- simulated_df_noreg(
                has_stock_init = TRUE,
                is_dominant = is_dominant,
                n_particle = n_particle,
                is_forecast = TRUE,
                cohort_type = cohort_type
            )
            do_forecast(
                df = df,
                work_dir = work_dir,
                quiet = TRUE
            )
            filename <- file.path(work_dir, "tmp-population-2000-Female.bin")
            con <- file(filename, "rb")
            n_meta <- readBin(con, what = "integer", n = 1)
            n_count <- readBin(con, what = "integer", n = 1)
            meta <- unserialize(readBin(con, what = "raw", n = n_meta))
            count <- readBin(con, what = "double", n = n_count)
            expect_identical(meta$cohort[[1]], 2000L)
            expect_identical(meta$sexgender[[1]], "Female")
            expect_false(anyNA(count))
            close(con)
            if (is_dominant) {
                filename <- file.path(work_dir, "tmp-births-2000-Female.bin")
                con <- file(filename, "rb")
                n_meta <- readBin(con, what = "integer", n = 1)
                n_count <- readBin(con, what = "integer", n = 1)
                meta <- unserialize(readBin(con, what = "raw", n = n_meta))
                count <- readBin(con, what = "double", n = n_count)
                expect_identical(meta$cohort[[1]], 2000L)
                expect_identical(meta$sexgender[[1]], "Female")
                expect_false(anyNA(count))
                close(con)
            }
            for (nm in c("deaths", "immigration1", "emigration1",
                         "immigration2", "emigration2")) {
                filename <- file.path(work_dir,
                                      sprintf("tmp-%s-2000-Female.bin", nm))
                con <- file(filename, "rb")
                n_meta <- readBin(con, what = "integer", n = 1)
                n_count <- readBin(con, what = "integer", n = 1)
                meta <- unserialize(readBin(con, what = "raw", n = n_meta))
                count <- readBin(con, what = "double", n = n_count)
                expect_identical(meta$cohort[[1]], 2000L)
                expect_identical(meta$sexgender[[1]], "Female")
                expect_false(anyNA(count))
                close(con)
            }
            diag <- readRDS(file.path(work_dir, "tmp-diagnostics-2000-Female.rds"))
            expect_false(anyNA(diag))
            unlink(work_dir, recursive = TRUE)
        }
    }
})
