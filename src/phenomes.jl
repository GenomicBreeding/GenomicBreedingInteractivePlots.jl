function sparsity(df_phenomes::DataFrame)::Matrix{Float64}
    .!Matrix(
        .!ismissing.(df_phenomes[:, 4:end]) .&& 
        .!isnan.(df_phenomes[:, 4:end]) .&& 
        .!isinf.(df_phenomes[:, 4:end]),
    )
end

function removesparsestroworcol(df_phenomes::DataFrame; prioritise_entries::Bool = true)::DataFrame
    S = sparsity(df_phenomes)
    # Remove sparsest row (entry) or column (trait)
    df_out = if prioritise_entries
        sparsity_rows = mean(S, dims = 2)[:, 1]
        if maximum(sparsity_rows) == 0.0
            df_phenomes
        else
            df_phenomes[sparsity_rows.<maximum(sparsity_rows), :]
        end
    else
        sparsity_cols = mean(S, dims = 1)[1, :]
        if maximum(sparsity_cols) == 0.0
            df_phenomes
        else
            df_phenomes[:, sparsity_cols.<maximum(sparsity_cols)]
        end
    end
    df_out
end

function loadphenomesdata(phenomes::Phenomes;
    idx_entries::Union{Nothing,Vector{Int64}} = nothing,
    idx_traits::Union{Nothing,Vector{Int64}} = nothing,
    threshold_n::Union{Nothing,Int64} = nothing,
    threshold_t::Union{Nothing,Int64} = nothing,
)::Tuple{DataFrame, Vector{String}, Vector{String}, Int64, Int64}
    # phenomes = Phenomes(n = 100, t = 3)
    # phenomes.entries = string.("entry_", 1:100)
    # phenomes.populations = StatsBase.sample(string.("pop_", 1:5), 100, replace = true)
    # phenomes.traits = ["trait_1", "trait_2", "long_trait_name number 3"]
    # phenomes.phenotypes = rand(Distributions.MvNormal([1, 2, 3], LinearAlgebra.I), 100)'
    # phenomes.phenotypes[1, 1] = missing
    # idx_entries = nothing; idx_traits = nothing; threshold_n = nothing; threshold_t = nothing
    # Check arguments
    if !checkdims(phenomes)
        throw(ArgumentError("The phenomes struct is corrupted."))
    end
    idx_entries = if isnothing(idx_entries)
        collect(1:length(phenomes.entries))
    else
        if (minimum(idx_entries) < 1) .|| maximum(idx_entries) > length(phenomes.entries)
            throw(
                ArgumentError(
                    "The indexes of the entries, `idx_entries` are out of bounds. Expected range: from 1 to " *
                    string(length(phenomes.entries)) *
                    " while the supplied range is from " *
                    string(minimum(idx_entries)) *
                    " to " *
                    string(maximum(idx_entries)) *
                    ".",
                ),
            )
        end
        idx_entries
    end
    idx_traits = if isnothing(idx_traits)
        collect(1:length(phenomes.traits))
    else
        if (minimum(idx_traits) < 1) .|| maximum(idx_traits) > length(phenomes.traits)
            throw(
                ArgumentError(
                    "The indexes of the traits, `idx_traits` are out of bounds. Expected range: from 1 to " *
                    string(length(phenomes.traits)) *
                    " while the supplied range is from " *
                    string(minimum(idx_traits)) *
                    " to " *
                    string(maximum(idx_traits)) *
                    ".",
                ),
            )
        end
        idx_traits
    end
    # Extract phenotypes and tabularise
    phenomes = if (length(idx_entries) < length(phenomes.entries)) && (length(idx_traits) < length(phenomes.traits))
        slice(phenomes, idx_entries = idx_entries, idx_traits = idx_traits)
    else
        phenomes
    end
    # Define the minimum number of entries and traits to retain
    threshold_n = if isnothing(threshold_n)
        # Default is 50% of the entries
        Int64(floor(length(phenomes.entries)*0.5))
    else
        if (threshold_n < 1) || (threshold_n > length(phenomes.entries))
            throw(
                ArgumentError(
                    "The threshold for the number of entries, `threshold_n` is out of bounds. Expected range: from 1 to " *
                    string(length(phenomes.entries)) *
                    " while the supplied value is " *
                    string(threshold_n) *
                    ".",
                ),
            )
        end
        threshold_n
    end
    threshold_t = if isnothing(threshold_t)
        # Default is 50% of the traits or 1
        length(phenomes.traits) == 1 ? 1 : Int64(floor(length(phenomes.traits)*0.5))
    else
        if (threshold_t < 1) || (threshold_t > length(phenomes.traits))
            throw(
                ArgumentError(
                    "The threshold for the number of traits, `threshold_t` is out of bounds. Expected range: from 1 to " *
                    string(length(phenomes.traits)) *
                    " while the supplied value is " *
                    string(threshold_t) *
                    ".",
                ),
            )
        end
        threshold_t
    end
    # Tabularise the phenomes
    df = tabularise(phenomes)
    traits = sort(unique(phenomes.traits))
    populations = sort(unique(phenomes.populations))
    for population in populations
        if sum(df.populations .== population) < 2
            throw(ArgumentError("Population $population has less than 2 entries."))
        end
    end
    # Output
    (
        df,
        traits,
        populations,
        threshold_n,
        threshold_t,
    )

end

function addpc1pc2(df_phenomes::DataFrame)::DataFrame
    traits = names(df)[4:end]
    t = length(traits)
    if t <= 2
        @warn "No need to do PCA as you only have at most 2 traits."
        return df
    end
    A::Matrix{Float64} = Matrix(df[:, 4:end])
    A = (A .- mean(A, dims = 1)) ./ std(A, dims = 1)
    # Remove traits with no variation
    v = var(A, dims = 1)[1, :]
    idx_cols = findall((abs.(v .- 1) .< 0.00001) .&& .!isnan.(v) .&& .!ismissing.(v) .&& .!isinf.(v))
    A = A[:, idx_cols]
    M = fit(PCA, A; maxoutdim = 2)
    vec_missing::Vector{Union{Missing,Float64}} = repeat([missing], size(df, 1))
    df.pc1 = deepcopy(vec_missing)
    df.pc2 = deepcopy(vec_missing)
    df[!, :pc1] = M.proj[:, 1]
    df[!, :pc2] = M.proj[:, 2]
    df
end

function plotinteractive2d(
    phenomes::Phenomes;
    idx_entries::Union{Nothing,Vector{Int64}} = nothing,
    idx_traits::Union{Nothing,Vector{Int64}} = nothing,
    threshold_n::Union{Nothing,Int64} = nothing,
    threshold_t::Union{Nothing,Int64} = nothing,
    prioritise_entries::Bool = true,
    impute::Bool = false,
)::Figure
    # phenomes = Phenomes(n = 100, t = 3)
    # phenomes.entries = string.("entry_", 1:100)
    # phenomes.populations = StatsBase.sample(string.("pop_", 1:5), 100, replace = true)
    # phenomes.traits = ["trait_1", "trait_2", "long_trait_name number 3"]
    # phenomes.phenotypes = rand(Distributions.MvNormal([1, 2, 3], LinearAlgebra.I), 100)'
    # n_missing = 20
    # phenomes.phenotypes[sample(1:length(phenomes.entries), n_missing, replace=true), sample(1:length(phenomes.traits), n_missing, replace=true)] .= missing
    # idx_entries = nothing; idx_traits = nothing; threshold_n = nothing; threshold_t = nothing; prioritise_entries = true; impute = true
    # Load phenomes data, tabularise, extract trait and population identifiers, minimum entry and trait thresholds while checking the arguments
    df, traits, populations, threshold_n, threshold_t = loadphenomesdata(phenomes, idx_entries = idx_entries, idx_traits = idx_traits)
    df = if !impute
        # Remove sparsest rows and columns until n and t become constant which means either:
        #   - no missing data remains
        #   - no data remains
        n = nrow(df)
        t = ncol(df) - 3
        for i in 1:maximum([length(phenomes.entries), length(phenomes.traits)])
            if (i > 1) && (nrow(df) == n) && ((ncol(df) - 3) == t)
                break
            end
            n -=  n - nrow(df)
            t -=  t - (ncol(df) - 3)
            df = removesparsestroworcol(df, prioritise_entries = prioritise_entries)
        end
        if n == 0
            throw(ArgumentError("Data is too sparse. No entries have been retained. Consider setting `impute` to true."))
        end
        if t == 0
            throw(ArgumentError("Data is too sparse. No traits have been retained. Consider setting `impute` to true."))
        end
        if n < threshold_n
            throw(ArgumentError("Data is too sparse. Number of entries is less than threshold_n ($n < $threshold_n). Consider setting `impute` to true."))
        end
        if t < threshold_t
            throw(ArgumentError("Data is too sparse. Number of traits is less than threshold_t ($t < $threshold_t). Consider setting `impute` to true."))
        end
        df
    else
        # Remove sparsest rows and columns until:
        #   - no missing data remains
        #   - number of entries is at least threshold_n
        #   - number of traits is at least threshold_t
        n = nrow(df)
        t = ncol(df) - 3
        for i in 1:maximum([length(phenomes.entries), length(phenomes.traits)])
            df = if (i > 1) && ((nrow(df) < n) || ((ncol(df) - 3) < t))
                break
            else
                removesparsestroworcol(df, prioritise_entries = prioritise_entries)
            end
            n -=  n - nrow(df)
            t -=  t - (ncol(df) - 3)
        end
        # Impute
        @warn "Imputing missing data via mean value imputation per trait."
        S = sparsity(df)
        sparsity_rows = mean(S, dims = 2)[:, 1]
        sparsity_cols = mean(S, dims = 1)[1, :]
        UnicodePlots.histogram(sparsity_rows, title = "Sparsity per entry")
        UnicodePlots.histogram(sparsity_cols, title = "Sparsity per locus-allele")
        for j in 4:ncol(df)
            # j = 4
            y = df[:, j]
            idx = .!ismissing.(y) .&& .!isnan.(y) .&& .!isinf.(y)
            μ_y = mean(y[idx])
            df[.!idx, j] .= μ_y
        end
        df
    end
    # Extract the first 2 principal components
    df = addpc1pc2(df)
    # Include the 2 PCs in the list of traits
    traits = names(df)[4:end]
    traits = traits[sortperm(uppercase.(traits))]
    # Activate Makie using the OpenGL plotting backend
    GLMakie.activate!()
    GLMakie.closeall() # close any open screen
    # Set plot size
    height, width = 1_200, 800
    # Define the entire plot/figure
    fig = begin
        fig = Figure(size = (height, width))
        # Set the trait 1 selector
        default_trait_x = if length(traits) > 2
            "pc1"
        else
            traits[1]
        end
        menu_trait_x = Menu(fig, options = traits, default = default_trait_x)
        tb_x = Textbox(fig)
        # Set the trait 2 selector
        default_trait_y = if length(traits) > 2
            "pc2"
        else
            traits[2]
        end
        menu_trait_y = Menu(fig, options = traits, default = default_trait_y)
        tb_y = Textbox(fig)
        # Instantiate Pearson's correlation value per pair of traits
        ρ = Observable(string(round(cor(df[!, default_trait_x], df[!, default_trait_y]), digits = 4)))
        # Place these trait menus in a left sidebar
        fig[1, 1] = vgrid!(
            Label(fig, "Search trait 1", width = nothing),
            tb_x,
            menu_trait_x,
            Label(fig, "Search trait 2", width = nothing),
            tb_y,
            menu_trait_y,
            Label(fig, @lift("ρ = $($ρ)"));
            tellheight = false,
        )
        # Place the main scatter plot with histograms for the 2 traits
        fig_main = fig[1:2, 2] = GridLayout()
        plot_hist_x = Axis(fig_main[1, 1])
        plot_scatter = Axis(fig_main[2, 1], xlabel = default_trait_x, ylabel = default_trait_y)
        plot_hist_y = Axis(fig_main[2, 2])

        plot_heatmap = Axis(
            fig[2, 1],
            title = "Trait correlations",
            xticks = (1:length(traits), traits),
            yticks = (1:length(traits), reverse(traits)),
            xaxisposition = :top,
            xticklabelrotation = deg2rad(90),
        )

        C = cor(Matrix(df[!, traits]))
        C = C[:, reverse(collect(1:end))]
        GLMakie.heatmap!(
            plot_heatmap,
            1:length(traits),
            1:length(traits),
            C,
            colorrange = (-1, 1),
            inspector_label = (self, (i, j), p) ->
                string("ρ = ", round(C[i, j], digits = 4), "\n", traits[i], "\n", reverse(traits)[j]),
        )


        X = Observable{Any}(df[!, default_trait_x])
        Y = Observable{Any}(df[!, default_trait_y])

        for pop in populations
            idx = findall(df.populations .== pop)

            x = @lift($X[idx])
            y = @lift($Y[idx])

            GLMakie.hist!(plot_hist_x, x)
            GLMakie.hist!(plot_hist_y, y, direction = :x)
            GLMakie.scatter!(
                plot_scatter,
                x,
                y,
                label = string(pop, " (n=", length(idx), ")"),
                inspector_label = (self, i, p) -> string(df.entries[idx][i], "\n(", df.populations[idx][i], ")"),
            )
        end
        connect!(ρ, @lift(rpad(round(cor($X, $Y), digits = 4), 2 * maximum(length.(traits)), " ")))

        leg = Legend(fig_main[1, 2], plot_scatter)
        GLMakie.hidedecorations!(plot_hist_x, grid = false)
        GLMakie.hidedecorations!(plot_hist_y, grid = false)
        leg.tellheight = true


        # Reactivity stuff
        # Search-ish functionality-ish
        on(tb_x.stored_string) do s
            bool_matches = .!isnothing.(match.(Regex(s, "i"), traits)) # case-insensitive with the Regex flag "i"
            if sum(bool_matches) > 0
                traits_reordered = vcat(traits[bool_matches], traits[.!bool_matches])
                menu_trait_x.selection = traits_reordered[1]
                menu_trait_x.i_selected = 1
                menu_trait_x.options = traits_reordered
                menu_trait_x.is_open = true
            else
                menu_trait_x.options = traits
            end
        end
        on(tb_y.stored_string) do s
            bool_matches = .!isnothing.(match.(Regex(s, "i"), traits)) # case-insensitive with the Regex flag "i"
            if sum(bool_matches) > 0
                traits_reordered = vcat(traits[bool_matches], traits[.!bool_matches])
                menu_trait_y.selection = traits_reordered[1]
                menu_trait_y.i_selected = 1
                menu_trait_y.options = traits_reordered
                menu_trait_y.is_open = true
            else
                menu_trait_y.options = traits
            end
        end
        # Drop down menus
        on(menu_trait_x.selection) do s
            X[] = df[!, s]
            plot_scatter.xlabel = s
            autolimits!(plot_scatter)
            autolimits!(plot_hist_x)
            autolimits!(plot_hist_y)
        end

        on(menu_trait_y.selection) do s
            # trait_y[] = s
            Y[] = df[!, s]
            plot_scatter.ylabel = s
            autolimits!(plot_scatter)
            autolimits!(plot_hist_x)
            autolimits!(plot_hist_y)
        end

        DataInspector(fig)
        fig
    end
    # Output
    fig
end
