function removesparsestroworcol(
    df_phenomes::DataFrame;
    prioritise_entries::Bool = true,
)::DataFrame
    S = .!Matrix(.!ismissing.(df_phenomes[:, 4:end]) .&& .!isnan.(df_phenomes[:, 4:end]) .&& .!isinf.(df_phenomes[:, 4:end]))
    # Remove sparsest row (entry) or column (trait)
    df_out = if prioritise_entries
        sparsity_rows = mean(S, dims=2)[:, 1]
        if maximum(sparsity_rows) == 0.0
            df_phenomes
        else
            df_phenomes[sparsity_rows .< maximum(sparsity_rows), :]
        end
    else
        sparsity_cols = mean(S, dims=1)[1, :]
        if maximum(sparsity_cols) == 0.0
            df_phenomes
        else
            df_phenomes[:, sparsity_cols .< maximum(sparsity_cols)]
        end
    end
    df_out
end


function plotinteractive2d(
    phenomes::Phenomes;
    idx_entries::Union{Nothing,Vector{Int64}} = nothing,
    idx_traits::Union{Nothing,Vector{Int64}} = nothing,
    prioritise_entries::Bool = true,
)::Figure
    # phenomes = Phenomes(n = 100, t = 3)
    # phenomes.entries = string.("entry_", 1:100)
    # phenomes.populations = StatsBase.sample(string.("pop_", 1:5), 100, replace = true)
    # phenomes.traits = ["trait_1", "trait_2", "long_trait_name number 3"]
    # phenomes.phenotypes = rand(Distributions.MvNormal([1, 2, 3], LinearAlgebra.I), 100)'
    # phenomes.phenotypes[1, 1] = missing
    # idx_entries = nothing; idx_traits = nothing; prioritise_entries = true
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
    df = tabularise(phenomes)
    traits = sort(unique(phenomes.traits))
    populations = sort(unique(phenomes.populations))
    for population in populations
        if sum(df.populations .== population) < 2
            throw(ArgumentError("Population $population has less than 2 entries."))
        end
    end

    n::Int64 = nrow(df)
    threshold_n = Int(floor(0.5*n))

    while (nrow(df) < n) || (nrow(df) <= threshold_n)
        n = nrow(df)
        df = removesparsestroworcol(df, prioritise_entries = prioritise_entries)
    end



    # Remove entries with at least 1 missing/NaN/Inf trait or remove trait/s with missing/NaN/Inf to keep all traits
    idx_rows, idx_cols = if prioritise_entries
        idx_rows = findall(
            mean(Matrix(.!ismissing.(df[:, 4:end]) .&& .!isnan.(df[:, 4:end]) .&& .!isinf.(df[:, 4:end])), dims = 2)[
                :,
                1,
            ] .== 1,
        )
        idx_cols = collect(1:size(df, 2))
        if length(idx_rows) == 0
            throw(
                ArgumentError(
                    "All entries have at least 1 missing trait value. Please consider:" * 
                    "\n\t(1) setting `prioritise_entries` to `false`," *
                    "\n\t(2) slicing the trial to include at least one non-sparse harvest, or" *
                    "\n\t(3) imputing missing phenotypes.",
                ),
            )
        end
        idx_rows, idx_cols
    else
        idx_rows = collect(1:size(df, 1))
        idx_cols = findall(
            mean(Matrix(.!ismissing.(df[:, 4:end]) .&& .!isnan.(df[:, 4:end]) .&& .!isinf.(df[:, 4:end])), dims = 1)[
                1,
                :,
            ] .== 1,
        )
        if length(idx_cols) == 0
            throw(
                ArgumentError(
                    "All traits have at least 1 entry with missing data. Please consider:" * 
                    "\n\t(1) setting `prioritise_entries` to `true`," *
                    "\n\t(2) slicing the trial to include at least one non-sparse harvest, or" *
                    "\n\t(3) imputing missing phenotypes.",
                ),
            )
        end
        idx_rows, idx_cols
    end
    df = df[idx_rows, idx_cols]






    # Extract the first 2 principal components
    traits = names(df)[4:end]
    t = length(traits)
    if t <= 2
        throw(ArgumentError("No need to do PCA as you only have at most 2 traits."))
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
