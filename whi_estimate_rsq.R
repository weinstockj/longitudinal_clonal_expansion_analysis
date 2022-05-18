suppressPackageStartupMessages({
    library(dplyr)
    library(magrittr)
    library(glue)
    library(tidyr)
    library(ggplot2)
    library(ggrepel)
    library(cowplot)
    library(stringr)
    library(rsample)
    library(purrr)
})

whi_clones = vroom::vroom("whi_clones.csv") %>%
    select(
        sample_id = `Sample.ID`,
        age = Age,
        CHR,
        position = `Position (hg19)`, # there are 4 rows with NA position!
        REF,
        ALT,
        VAF,
        high_confidence = `High.Conf`
    ) %>%
    dplyr::filter(high_confidence == "yes") %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(
        variant_label = glue::glue("{CHR}-{position}-{REF}-{ALT}"),
        n_drivers = length(unique(variant_label))
    ) %>%
    dplyr::filter(n_drivers == 1) %>%
    dplyr::filter(n() >= 3) %>%
    dplyr::ungroup(.)


training = whi_clones %>%
    group_by(sample_id) %>%
    arrange(age) %>%
    drop_na %>%
    slice(1:(n() - 1))  %>%# exclude last timepoint
    ungroup %>%
    mutate(id = 1:n()) %>%
    group_by(sample_id) %>%
    tidyr::nest(.) %>%
    mutate(
       unpooled_fitness_model = map(
            data,
            ~lm(VAF ~ age, data = .x)
       ),
       unpooled_fitness_estimate = map_dbl(
            unpooled_fitness_model,
            ~broom::tidy(.x) %>% dplyr::filter(term == "age") %>% dplyr::pull(estimate)
       )
    ) %>%
    tidyr::unnest(data) %>%
    dplyr::ungroup(.) # 75 rows from 23 individuals

testing = whi_clones %>% # 23 rows
    group_by(sample_id) %>%
    arrange(age) %>%
    drop_na %>%
    slice(n()) %>% # extract last timepoint for testing
    ungroup %>%
    mutate(is_testing = TRUE)

fitness_estimates = training %>%
            select(sample_id, unpooled_fitness_estimate) %>%
            distinct(.)

dfm = testing %>%
    bind_rows(training %>% dplyr::select(-matches("unpooled"))) %>%
    group_by(sample_id) %>%
    arrange(age) %>%
    mutate(
        prev_VAF = lag(VAF),
        prev_age = lag(age),
        dVAF = VAF - prev_VAF,
        dT = age - prev_age,
        dVAFdT = dVAF / dT,
        dVAFdT = pmax(0, dVAFdT) # truncate at 0
    ) %>%
    ungroup %>%
    dplyr::inner_join(fitness_estimates) %>%
    tidyr::drop_na(is_testing) # filter to test data points

dfm_bt = dfm %>%
    bootstraps(times = 1000) %>%
    mutate(
        unpooled_fitness_rsq = map_dbl(
            splits,
            ~with(as.data.frame(.x), cor(unpooled_fitness_estimate, inverse_normalize(dVAFdT), use = "complete.obs") ^ 2)
        )
    )

mean_rsquared = dfm_bt %>%
    pull(unpooled_fitness_rsq) %>%
    mean
pacer_rsquared = .325

bt_plot = ggplot(data = dfm_bt, aes(x = unpooled_fitness_rsq)) +
            geom_histogram() +
            geom_vline(xintercept = pacer_rsquared, linetype = "dashed", color = "blue") +
            geom_vline(xintercept = mean_rsquared, linetype = "dashed", color = "red") +
            scale_x_continuous(breaks = scales::pretty_breaks(n = 8), labels = scales::percent) + 
            scale_y_continuous(expand = c(0, 0)) + 
            cowplot::theme_cowplot(font_size = 12) +
            labs(
                x = "Boostrapped Rsq estimates\nfrom unpooled linear model predictions of dVAFdT",
                y = "Count per 1,000 bootstrap replicates"
            )

date = "2022_05_13"
ggsave(glue("WHI_unpooled_linear_dVAFdT_hist_{date}.pdf"), bt_plot, width = 6, height = 4, units = "in")
