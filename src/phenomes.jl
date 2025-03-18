"""
    loadphenomesdata(phenomes::Phenomes; 
                     idx_entries::Union{Nothing,Vector{Int64}} = nothing,
                     idx_traits::Union{Nothing,Vector{Int64}} = nothing, 
                     threshold_n::Union{Nothing,Int64} = nothing,
                     threshold_t::Union{Nothing,Int64} = nothing
                    )::Tuple{DataFrame,Vector{String},Int64,Int64}

Load and process phenotypic data from a Phenomes struct into a tabular format.

# Arguments
- `phenomes::Phenomes`: A Phenomes struct containing phenotypic data
- `idx_entries::Union{Nothing,Vector{Int64}}=nothing`: Optional vector of indices to select specific entries
- `idx_traits::Union{Nothing,Vector{Int64}}=nothing`: Optional vector of indices to select specific traits
- `threshold_n::Union{Nothing,Int64}=nothing`: Minimum number of entries to retain (defaults to 50% of entries)
- `threshold_t::Union{Nothing,Int64}=nothing`: Minimum number of traits to retain (defaults to 50% of traits or 1)

# Returns
A tuple containing:
- `DataFrame`: Tabularized phenotypic data
- `Vector{String}`: Sorted unique trait names
- `Int64`: Final threshold value for number of entries
- `Int64`: Final threshold value for number of traits

# Throws
- `ArgumentError`: If phenomes struct is corrupted
- `ArgumentError`: If idx_entries are out of bounds
- `ArgumentError`: If idx_traits are out of bounds
- `ArgumentError`: If threshold_n is out of bounds
- `ArgumentError`: If threshold_t is out of bounds
- `ArgumentError`: If any population has less than 2 entries

# Example
```jldoctest; setup = :(using GBCore, GBPlotsInteractive, StatsBase, LinearAlgebra, Distributions)
julia> phenomes = Phenomes(n = 100, t = 3);

julia> phenomes.entries = string.("entry_", 1:100);

julia> phenomes.populations = StatsBase.sample(string.("pop_", 1:5), 100, replace = true);

julia> phenomes.traits = ["trait_1", "trait_2", "long_trait_name number 3"];

julia> phenomes.phenotypes = rand(Distributions.MvNormal([1, 2, 3], LinearAlgebra.I), 100)';

julia> n_missing = 20;

julia> phenomes.phenotypes[sample(1:length(phenomes.entries), n_missing, replace=true), sample(1:length(phenomes.traits), n_missing, replace=true)] .= missing;

julia> (df, traits, threshold_n, threshold_t) = loadphenomesdata(phenomes);

julia> size(df)
(100, 6)

julia> length(traits)
3

julia> (threshold_n, threshold_t)
(50, 1)
```
"""
function loadphenomesdata(
    phenomes::Phenomes;
    idx_entries::Union{Nothing,Vector{Int64}} = nothing,
    idx_traits::Union{Nothing,Vector{Int64}} = nothing,
    threshold_n::Union{Nothing,Int64} = nothing,
    threshold_t::Union{Nothing,Int64} = nothing,
)::Tuple{DataFrame,Vector{String},Int64,Int64}
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
        Int64(floor(length(phenomes.entries) * 0.5))
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
        length(phenomes.traits) == 1 ? 1 : Int64(floor(length(phenomes.traits) * 0.5))
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
    (df, traits, threshold_n, threshold_t)

end

"""
    sparsity(df::DataFrame)::Matrix{Float64}

Calculate the sparsity pattern of a DataFrame's numerical columns.

Takes a DataFrame and returns a boolean matrix where:
- `true` indicates missing, NaN, or infinite values
- `false` indicates valid numerical values

The function considers only columns from the 4th column onwards.

# Arguments
- `df::DataFrame`: Input DataFrame to analyze

# Returns
- `Matrix{Float64}`: A boolean matrix where `true` represents missing/invalid data points

# Examples
"""
function sparsity(df::DataFrame)::Matrix{Float64}
    .!Matrix(.!ismissing.(df[:, 4:end]) .&& .!isnan.(df[:, 4:end]) .&& .!isinf.(df[:, 4:end]))
end

"""
    removesparsestroworcol(df::DataFrame; prioritise_entries::Bool = true)::DataFrame

Remove either the sparsest row or column from a DataFrame based on sparsity values.

# Arguments
- `df::DataFrame`: Input DataFrame to process
- `prioritise_entries::Bool = true`: If true, prioritizes removing sparse columns (loci-alleles). If false, prioritizes removing sparse rows (entries).

# Returns
- `DataFrame`: A new DataFrame with the sparsest row or column removed. If no sparsity is found (all values are 0.0), returns the original DataFrame unchanged.

# Details
The function calculates sparsity for each row or column (depending on `prioritise_entries`), and removes the one with the highest sparsity value. Sparsity is computed as the mean of sparsity values across the respective dimension.

# Example
"""
function removesparsestroworcol(df::DataFrame; prioritise_entries::Bool = true)::DataFrame
    S = sparsity(df)
    # Remove sparsest row (entry) or column (trait)
    df_out = if !prioritise_entries
        sparsity_rows = mean(S, dims = 2)[:, 1]
        if maximum(sparsity_rows) == 0.0
            df
        else
            df[sparsity_rows.<maximum(sparsity_rows), :]
        end
    else
        sparsity_cols = mean(S, dims = 1)[1, :]
        idx_cols = vcat(collect(1:3), collect(4:ncol(df))[sparsity_cols.<maximum(sparsity_cols)])
        if maximum(sparsity_cols) == 0.0
            df
        else
            df[:, idx_cols]
        end
    end
    df_out
end

"""
    filterphenomesdata(df::DataFrame; 
        threshold_n::Union{Nothing,Int64} = nothing,
        threshold_t::Union{Nothing,Int64} = nothing,
        prioritise_entries::Bool = true,
        impute::Bool = false)::DataFrame

Filter and process phenome data from a DataFrame by handling missing values through either removal or imputation.

# Arguments
- `df::DataFrame`: Input DataFrame containing phenomes data with entries in rows and traits in columns (from column 4 onwards)
- `threshold_n::Union{Nothing,Int64}`: Minimum number of entries required in the filtered dataset
- `threshold_t::Union{Nothing,Int64}`: Minimum number of traits required in the filtered dataset
- `prioritise_entries::Bool`: If true, prioritizes keeping entries over traits when removing sparse data
- `impute::Bool`: If true, performs mean value imputation instead of removing sparse data

# Returns
- `DataFrame`: Filtered and processed DataFrame

# Details
When `impute=false`:
- Iteratively removes sparsest rows or columns until no missing data remains or no data is left
- Throws ArgumentError if resulting data doesn't meet threshold requirements

When `impute=true`:
- Removes sparsest rows and columns until meeting threshold requirements
- Performs mean value imputation for remaining missing values
- Displays histograms of sparsity distribution

# Throws
- `ArgumentError`: If resulting data is too sparse or doesn't meet threshold requirements when `impute=false`

# Note
The input DataFrame should have entry information in the first three columns and traits data from column 4 onwards.
"""
function filterphenomesdata(
    df::DataFrame;
    threshold_n::Union{Nothing,Int64} = nothing,
    threshold_t::Union{Nothing,Int64} = nothing,
    prioritise_entries::Bool = true,
    impute::Bool = false,
)::DataFrame
    entries = unique(df.entries)
    traits = names(df)[4:end]
    df = if !impute
        # Remove sparsest rows and columns until n and t become constant which means either:
        #   - no missing data remains
        #   - no data remains
        n = length(entries)
        t = length(traits)
        for i = 1:maximum([length(entries), length(traits)])
            if (i > 1) && (nrow(df) == n) && ((ncol(df) - 3) == t)
                break
            end
            n -= n - nrow(df)
            t -= t - (ncol(df) - 3)
            df = removesparsestroworcol(df, prioritise_entries = prioritise_entries)
        end
        if n == 0
            throw(ArgumentError("Data is too sparse. No entries have been retained. Consider setting `impute` to true."))
        end
        if t == 0
            throw(ArgumentError("Data is too sparse. No traits have been retained. Consider setting `impute` to true."))
        end
        if n < threshold_n
            throw(
                ArgumentError(
                    "Data is too sparse. Number of entries is less than threshold_n ($n < $threshold_n). Consider setting `impute` to true.",
                ),
            )
        end
        if t < threshold_t
            throw(
                ArgumentError(
                    "Data is too sparse. Number of traits is less than threshold_t ($t < $threshold_t). Consider setting `impute` to true.",
                ),
            )
        end
        df
    else
        # Remove sparsest rows and columns until:
        #   - no missing data remains
        #   - number of entries is at least threshold_n
        #   - number of traits is at least threshold_t
        n = nrow(df)
        t = ncol(df) - 3
        for i = 1:maximum([length(entries), length(traits)])
            # i = 1
            df = if (i > 1) && ((nrow(df) < n) || ((ncol(df) - 3) < t))
                break
            else
                removesparsestroworcol(df, prioritise_entries = prioritise_entries)
            end
            n -= n - nrow(df)
            t -= t - (ncol(df) - 3)
        end
        # Impute
        @warn "Imputing missing data via mean value imputation per trait."
        S = sparsity(df)
        sparsity_rows = mean(S, dims = 2)[:, 1]
        sparsity_cols = mean(S, dims = 1)[1, :]
        UnicodePlots.histogram(sparsity_rows, title = "Sparsity per entry")
        UnicodePlots.histogram(sparsity_cols, title = "Sparsity per locus-allele")
        for j = 4:ncol(df)
            # j = 4
            y = df[:, j]
            idx = .!ismissing.(y) .&& .!isnan.(y) .&& .!isinf.(y)
            μ_y = mean(y[idx])
            df[.!idx, j] .= μ_y
        end
        df
    end
    # Output
    df
end

"""
    addpc1pc2(df::DataFrame)::DataFrame

Compute and append the first two principal components (PC1 and PC2) of the phenotype data

This function performs Principal Component Analysis (PCA) on the trait columns of the input DataFrame
(assumed to start from column 4 onwards). The traits are standardized before PCA is performed.

# Arguments
- `df::DataFrame`: Input DataFrame with trait columns starting from column 4

# Returns
- `DataFrame`: Original DataFrame with two additional columns:
  - `pc1`: First principal component scores
  - `pc2`: Second principal component scores

# Notes
- If there are 2 or fewer traits, a warning is issued and the original DataFrame is returned unchanged
- Traits with no variation are removed before PCA
- The trait data is standardized (centered and scaled) before PCA
- Missing, infinite, or NaN values in variance are handled

# Examples
"""
function addpc1pc2(df::DataFrame)::DataFrame
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

"""
    figurelayout(df::DataFrame; height::Int64=1_200, width::Int64=800, traits::Vector{String})
    -> Tuple{Figure,GridLayout,Axis,Axis,Axis,Axis,Observable{String},String,String,Menu,Menu,Textbox,Textbox}

Create an interactive figure layout for visualizing trait relationships in a DataFrame.

# Arguments
- `df::DataFrame`: Input DataFrame containing trait data
- `height::Int64=1_200`: Height of the figure in pixels
- `width::Int64=800`: Width of the figure in pixels  
- `traits::Vector{String}`: Vector of trait names to be visualized

# Returns
Tuple containing:
- `Figure`: Main figure object
- `GridLayout`: Layout for the main scatter plot and histograms
- `Axis`: Histogram axis for trait 1
- `Axis`: Scatter plot axis
- `Axis`: Histogram axis for trait 2  
- `Axis`: Heatmap axis for trait correlations
- `Observable{String}`: Observable containing Pearson correlation value
- `String`: Default trait 1 name
- `String`: Default trait 2 name
- `Menu`: Dropdown menu for selecting trait 1
- `Menu`: Dropdown menu for selecting trait 2
- `Textbox`: Search box for trait 1
- `Textbox`: Search box for trait 2

The function creates an interactive visualization layout with:
- A scatter plot showing relationship between two selected traits
- Marginal histograms for each trait
- A heatmap showing correlations between all traits
- Interactive menus to select traits
- Search boxes to find traits
- Display of Pearson correlation between selected traits
"""
function figurelayout(
    df::DataFrame;
    height::Int64 = 1_200,
    width::Int64 = 800,
    traits::Vector{String},
)::Tuple{Figure,GridLayout,Axis,Axis,Axis,Axis,Observable{String},String,String,Menu,Menu,Textbox,Textbox}
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
    # Set-up the heatmap of trait correlations
    plot_heatmap = Axis(
        fig[2, 1],
        title = "Trait correlations",
        xticks = (1:length(traits), traits),
        yticks = (1:length(traits), reverse(traits)),
        xaxisposition = :top,
        xticklabelrotation = deg2rad(90),
    )
    # Output
    (
        fig,
        fig_main,
        plot_hist_x,
        plot_scatter,
        plot_hist_y,
        plot_heatmap,
        ρ,
        default_trait_x,
        default_trait_y,
        menu_trait_x,
        menu_trait_y,
        tb_x,
        tb_y,
    )
end

"""
    heatmapinteractive!(plot_heatmap::Axis; df::DataFrame, traits::Vector{String})

Creates an interactive correlation heatmap for the specified traits using GLMakie.

# Arguments
- `plot_heatmap::Axis`: The axis object where the heatmap will be plotted
- `df::DataFrame`: DataFrame containing the data for correlation analysis
- `traits::Vector{String}`: Vector of trait names to be included in the correlation analysis

# Details
- Calculates correlation matrix for specified traits from the DataFrame
- Reverses column order for display purposes
- Creates interactive heatmap with correlation values ranging from -1 to 1
- Includes inspector label showing correlation coefficient and trait pairs on hover

# Returns
Nothing, modifies the plot in-place.
"""
function heatmapinteractive!(plot_heatmap::Axis; df::DataFrame, traits::Vector{String})
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
end

"""
    scatterplotinteractive!(X::Observable, Y::Observable, ρ::Observable, plot_scatter::Axis, 
                           plot_hist_x::Axis, plot_hist_y::Axis; df::DataFrame, traits::Vector{String})

Create an interactive scatter plot with marginal histograms for multiple populations.

# Arguments
- `X::Observable`: Observable containing x-axis (trait 1) data
- `Y::Observable`: Observable containing y-axis (trait 2) data
- `ρ::Observable`: Observable for Pearson's correlation coefficient for `x` and `y` traits to display
- `plot_scatter::Axis`: Main scatter plot axis
- `plot_hist_x::Axis`: X-axis marginal histogram axis
- `plot_hist_y::Axis`: Y-axis marginal histogram axis
- `df::DataFrame`: DataFrame containing the data with required columns:
    - `populations`: Population identifiers
    - `entries`: Entry labels
- `traits::Vector{String}`: Vector of trait names for formatting

# Details
Creates a scatter plot with population-specific colours and interactive labels. Includes marginal
histograms for both x and y dimensions. Updates correlation coefficient dynamically.

# Returns
Nothing, modifies the input axes in place.
"""
function scatterplotinteractive!(
    X::Observable,
    Y::Observable,
    ρ::Observable,
    plot_scatter::Axis,
    plot_hist_x::Axis,
    plot_hist_y::Axis;
    df::DataFrame,
    traits::Vector{String},
)
    populations = sort(unique(df.populations))
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
    GLMakie.hidedecorations!(plot_hist_x, grid = false)
    GLMakie.hidedecorations!(plot_hist_y, grid = false)
    connect!(ρ, @lift(rpad(round(cor($X, $Y), digits = 4), 2 * maximum(length.(traits)), " ")))
end

"""
    plotinteractive2d(phenomes::Phenomes; 
        idx_entries::Union{Nothing,Vector{Int64}} = nothing,
        idx_traits::Union{Nothing,Vector{Int64}} = nothing,
        threshold_n::Union{Nothing,Int64} = nothing,
        threshold_t::Union{Nothing,Int64} = nothing,
        prioritise_entries::Bool = true,
        impute::Bool = false
    )::Figure

Creates an interactive 2D visualization of phenotypic data with the following components:
- A scatter plot of two selected traits
- Marginal histograms for both traits
- A correlation heatmap of all traits

# Arguments
- `phenomes::Phenomes`: A Phenomes object containing the phenotypic data
- `idx_entries::Union{Nothing,Vector{Int64}}`: Optional indices to subset entries
- `idx_traits::Union{Nothing,Vector{Int64}}`: Optional indices to subset traits
- `threshold_n::Union{Nothing,Int64}`: Minimum number of non-missing values per entry
- `threshold_t::Union{Nothing,Int64}`: Minimum number of non-missing values per trait
- `prioritise_entries::Bool`: If true, prioritizes keeping entries over traits when filtering
- `impute::Bool`: If true, performs imputation on missing values

# Returns
- `Figure`: A Makie figure object containing the interactive visualization

# Features
- Interactive trait selection via dropdown menus
- Search functionality for traits
- Data inspection on hover
- Automatic correlation calculation
- Principal component analysis visualization
- Responsive layout with linked histograms and scatter plot

# Notes
- Uses GLMakie as the backend
- Automatically handles missing data through filtering and optional imputation
- Includes first two principal components in the trait selection options
"""
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
    df, traits, threshold_n, threshold_t = loadphenomesdata(
        phenomes,
        idx_entries = idx_entries,
        idx_traits = idx_traits,
        threshold_n = threshold_n,
        threshold_t = threshold_t,
    )
    # Filter phenomes data until no missing data remains or until the minimum thresholds are reached and then imputed if necessary
    df = filterphenomesdata(
        df,
        threshold_n = threshold_n,
        threshold_t = threshold_t,
        prioritise_entries = prioritise_entries,
        impute = impute,
    )
    # Extract the first 2 principal components
    df = addpc1pc2(df)
    # Include the 2 PCs in the list of traits and sort alphabetically
    traits = names(df)
    traits = traits[isnothing.(match.(Regex("id|entries|populations", "i"), traits))]
    traits = traits[sortperm(uppercase.(traits))]
    # Activate Makie using the OpenGL plotting backend and close any open screen
    GLMakie.activate!()
    GLMakie.closeall()
    # Set-up the figure layout
    (
        fig,
        fig_main,
        plot_hist_x,
        plot_scatter,
        plot_hist_y,
        plot_heatmap,
        ρ,
        default_trait_x,
        default_trait_y,
        menu_trait_x,
        menu_trait_y,
        tb_x,
        tb_y,
    ) = figurelayout(df, traits = traits)
    # Compute the correlation matrix and plot the heatmap
    heatmapinteractive!(plot_heatmap, df = df, traits = traits)
    # Plot the scatter plot with histograms
    X = Observable{Any}(df[!, default_trait_x])
    Y = Observable{Any}(df[!, default_trait_y])
    scatterplotinteractive!(X, Y, ρ, plot_scatter, plot_hist_x, plot_hist_y, df = df, traits = traits)
    # Add the legend
    leg = Legend(fig_main[1, 2], plot_scatter)
    leg.tellheight = true
    # Reactivity stuff
    ## Hacky search functionality
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
    ## Drop down menus
    on(menu_trait_x.selection) do s
        X[] = df[!, s]
        plot_scatter.xlabel = s
        autolimits!(plot_scatter)
        autolimits!(plot_hist_x)
        autolimits!(plot_hist_y)
    end
    on(menu_trait_y.selection) do s
        Y[] = df[!, s]
        plot_scatter.ylabel = s
        autolimits!(plot_scatter)
        autolimits!(plot_hist_x)
        autolimits!(plot_hist_y)
    end
    # Add a data inspector
    DataInspector(fig)
    # Plot
    fig
end
