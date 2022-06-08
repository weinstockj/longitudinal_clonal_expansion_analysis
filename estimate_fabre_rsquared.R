# library(brms)
library(dplyr)
library(logger)
library(stringr)
library(glue)
library(tidyr)
library(rsample)
library(purrr)
library(rstan)
library(tidybayes)
# Uses SardID, and needed to get depths
raw = vroom::vroom("clonal_dynamics/data/ALLvariants_exclSynonymous_Xadj.txt")  

# Uses PID
s2 = readxl::read_excel("media-12.xlsx", skip = 3) 
# their 15 validation samples, not used at the moment since mostly splicing
s4 = readxl::read_excel("media-14.xlsx", skip = 3) %>%
            rename(
                PID = `Sample ID`,
                protein = Protein
            ) %>%
            mutate(protein = str_remove(protein, "p."))

mergeds2 = raw %>% 
    select(Sex = Gender, `Study Phase` = Phase, Age, Gene, Start = START, End = END, SardID, MUTcount = MUTcount_Xadj, depth = TOTALcount) %>%  
    inner_join(
          s2 
    ) %>%
    rename(PID = `Sample ID`, age = Age, protein = Protein) %>%
    mutate(
        protein = str_remove(protein, "p."),
        gene_protein = glue("{Gene}-{protein}"), 
        age = as.integer(age) # needed to merge on this column
    )

annotated = vroom::vroom("sidd_fabre_clones.csv") %>%
            select(
                PID = `Sample ID`,
                Gene,
                protein = Protein,
                age = Age,
                High.Conf
            ) %>%
            mutate(
                protein = str_remove(protein, "p."),
                High.Conf = High.Conf == "yes",
                gene_protein = glue("{Gene}-{protein}"),
                age = as.integer(age) #  converting to integer helps to join dataframes
            ) %>%
            left_join( # merged in depths
                mergeds2 %>%
                    select(PID, age, gene_protein, MUTcount, depth, VAF)
            ) %>%
            filter(High.Conf) %>%
            add_count(PID, gene_protein, name = "n_timepoints") %>%
            group_by(PID) %>%
            mutate(n_drivers = length(unique(gene_protein))) %>%
            ungroup %>%
            filter(n_timepoints >= 3 & n_drivers == 1) # 43 samples have 2 or fewer timepoints and 2141 and 179 out of 285 donors have 2 or more drivers

training = annotated %>%
    group_by(PID) %>%
    arrange(age) %>%
    drop_na %>%
    slice(1:(n() - 1))  %>%# exclude last timepoint
    ungroup %>%
    mutate(id = 1:n()) %>%
    group_by(PID) %>%
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
    dplyr::ungroup(.)

testing = annotated %>%
    group_by(PID) %>%
    arrange(age) %>%
    drop_na %>%
    slice(n()) %>% # extract last timepoint for testing
    ungroup


#https://cran.r-project.org/web/packages/brms/vignettes/brms_customfamilies.html#the-beta-binomial-distribution
#
# stan_funs <- "
#     real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
#         return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
#     }
#     int beta_binomial2_rng(real mu, real phi, int T) {
#         return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
#     }
# "
# stanvars <- stanvar(scode = stan_funs, block = "functions") 
#
# prior = c(
#     prior(cauchy(0, 0.05), class = sd, coef = age, group = PID),
#     prior(cauchy(0, 0.10), class = sd, coef = age, group = Gene),
#     prior(cauchy(0, 0.10), class = sd, coef = age, group = protein)
#     # prior(uniform(-30, 0), class = b, coef = Intercept, group = id)
# )
#



log_info("fitting model now")

# m = stan(
#     model_code = readLines("stancode.txt"),
#     chains = 5,
#     iter = 800,
#     cores = 5,
#     data = tibble::lst(
#        N = nrow(training),
#        Y = as.integer(training$MUTcount),
#        depth = as.integer(training$depth),
#        beta = 218L,
#        age = training$age,
#        N_1 = length(unique(training$Gene)),
#        N_2 = length(unique(training$protein)),
#        N_3 = length(unique(training$PID)),
#        N_4 = N_3,
#        J_1 = as.integer(as.factor(training$Gene)),
#        J_2 = as.integer(as.factor(training$protein)),
#        J_3 = as.integer(as.factor(training$PID)),
#        J_4 = J_3,
#        prior_only = FALSE,
#
#        N_new = nrow(testing),
#        VAF_new = testing$VAF,
#        depth_new = as.integer(testing$depth),
#        age_new = testing$age,
#        J_1_new = as.integer(factor(testing$Gene, levels = levels(as.factor(training$Gene)))),
#        J_2_new = as.integer(factor(testing$protein, levels = levels(as.factor(training$protein)))),
#        J_3_new = as.integer(factor(testing$PID, levels = levels(as.factor(training$PID)))),
#        J_4_new = J_3_new
#     )
# )
log_info("done")

lookup = tibble(
    # .variable = glue("mu_new[{1:nrow(testing)}]"),
    .variable = glue("total_effect[{1:nrow(testing)}]"),
    row = 1:nrow(testing)
)

param_lookup = bind_rows(
    tibble(
        label = unique(training$Gene)
    ) %>%
        mutate(index = 1:n(), .variable = glue("gene[{index}]")),
    tibble(
        label = unique(training$protein)
    ) %>%
        mutate(index = 1:n(), .variable = glue("site[{index}]")),
    tibble(
        label = unique(training$PID)
    ) %>%
        mutate(index = 1:n(), .variable = glue("PID[{index}]"))
)

inverse_normalize = function(x) {
    qnorm(rank(x, na.last = "keep") / (sum(!is.na(x)) + 1))
}

clone_fitness = tidy_draws(m) %>%
    gather_draws(`total_effect.*`, regex = TRUE) %>%
    mean_hdi  %>%
    distinct(.variable, .value) %>%
    inner_join(lookup) %>%
    inner_join(testing %>% mutate(row = 1:nrow(testing))) %>%
    rename(fabre_fitness = `.value`) 

annotated_with_fitness = clone_fitness %>%
    bind_rows(training %>% dplyr::select(-matches("unpooled"))) %>%
    group_by(PID) %>%
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
    dplyr::left_join(
        training %>% 
            dplyr::select(PID, unpooled_fitness_estimate) %>%
            dplyr::distinct(.)
    )

readr::write_tsv(annotated_with_fitness, "fabre_with_total_fitness_estimate.tsv")
model = lm(inverse_normalize(dVAFdT) ~ fabre_fitness, data = annotated_with_fitness)
summary(model) # Rsq of .03473

dfm_bt = annotated_with_fitness %>%
    tidyr::drop_na(fabre_fitness) %>% # filter to test data points
    bootstraps(times = 1000) %>%
    mutate(
        # models = map(splits, ~lm(inverse_normalize(dVAFdT) ~ dVAFdT_pred, data = .x)),
        # tidy = map(models, ~broom::glance(.x)),
        models_fitness = map(splits, ~lm(inverse_normalize(dVAFdT) ~ fabre_fitness, data = .x)),
        tidy_fitness = map(models_fitness, ~broom::glance(.x)),
        unpooled_fitness_rsq = map_dbl(
            splits,
            ~with(as.data.frame(.x), cor(unpooled_fitness_estimate, inverse_normalize(dVAFdT), use = "complete.obs") ^ 2)
        )
    )

fabre_mean_rsquared = dfm_bt %>%
    unnest(tidy_fitness) %>%
    pull(r.squared) %>%
    mean
pacer_rsquared = .325

bt_plot = ggplot(data = dfm_bt %>% unnest(tidy_fitness), aes(x = r.squared)) +
            geom_histogram() +
            geom_vline(xintercept = pacer_rsquared, linetype = "dashed", color = "blue") +
            geom_vline(xintercept = fabre_mean_rsquared, linetype = "dashed", color = "red") +
            scale_x_continuous(breaks = scales::pretty_breaks(n = 8), labels = scales::percent) + 
            scale_y_continuous(expand = c(0, 0)) + 
            cowplot::theme_cowplot(font_size = 12) +
            labs(
                x = "Boostrapped Rsq estimates\nfrom Fabre et al. predictions of dVAFdT",
                y = "Count per 1,000 bootstrap replicates"
            )

date = "2022_06_08"
ggsave(glue("expanded_fabre_dVAFdT_hist_{date}.pdf"), bt_plot, width = 6, height = 4, units = "in")

unpooled_mean_rsquared = dfm_bt %>%
    pull(unpooled_fitness_rsq) %>%
    mean

bt_unpooled_plot = ggplot(data = dfm_bt, aes(x = unpooled_fitness_rsq)) +
            geom_histogram() +
            geom_vline(xintercept = pacer_rsquared, linetype = "dashed", color = "blue") +
            geom_vline(xintercept = unpooled_mean_rsquared, linetype = "dashed", color = "red") +
            scale_x_continuous(breaks = scales::pretty_breaks(n = 8), labels = scales::percent) + 
            scale_y_continuous(expand = c(0, 0)) + 
            cowplot::theme_cowplot(font_size = 12) +
            labs(
                x = "Boostrapped Rsq estimates\nfrom unpooled linear model of dVAFdT",
                y = "Count per 1,000 bootstrap replicates"
            )

ggsave(glue("expanded_unpooled_linear_dVAFdT_hist_{date}.pdf"), bt_unpooled_plot, width = 6, height = 4, units = "in")
