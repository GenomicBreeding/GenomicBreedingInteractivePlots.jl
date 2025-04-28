function loadgenomesdata(
    genomes::Phenomes;
    idx_entries::Union{Nothing,Vector{Int64}} = nothing,
    idx_loci_alleles::Union{Nothing,Vector{Int64}} = nothing,
    threshold_n::Union{Nothing,Int64} = nothing,
    threshold_p::Union{Nothing,Int64} = nothing,
    n_pcs::Int64 = 10,
)
    # )::Tuple{DataFrame,Vector{String},Int64,Int64}
    # genomes = simulategenomes(sparsity=0.1, verbose=false)
    # idx_entries = nothing; idx_loci_alleles = nothing; threshold_n = nothing; threshold_p = nothing
    # Check arguments
    if !checkdims(genomes)
        throw(ArgumentError("The genomes struct is corrupted ☹."))
    end
    idx_entries = if isnothing(idx_entries)
        collect(1:length(genomes.entries))
    else
        if (minimum(idx_entries) < 1) .|| maximum(idx_entries) > length(genomes.entries)
            throw(
                ArgumentError(
                    "The indexes of the entries, `idx_entries` are out of bounds. Expected range: from 1 to " *
                    string(length(genomes.entries)) *
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
    idx_loci_alleles = if isnothing(idx_loci_alleles)
        collect(1:length(genomes.loci_alleles))
    else
        if (minimum(idx_loci_alleles) < 1) .|| maximum(idx_loci_alleles) > length(genomes.loci_alleles)
            throw(
                ArgumentError(
                    "The indexes of the traits, `idx_loci_alleles` are out of bounds. Expected range: from 1 to " *
                    string(length(genomes.loci_alleles)) *
                    " while the supplied range is from " *
                    string(minimum(idx_loci_alleles)) *
                    " to " *
                    string(maximum(idx_loci_alleles)) *
                    ".",
                ),
            )
        end
        idx_loci_alleles
    end
    # Slice genomes
    genomes =
        if (length(idx_entries) < length(genomes.entries)) && (length(idx_loci_alleles) < length(genomes.loci_alleles))
            slice(genomes, idx_entries = idx_entries, idx_loci_alleles = idx_loci_alleles)
        else
            genomes
        end
    # Define the minimum number of entries and traits to retain
    threshold_n = if isnothing(threshold_n)
        # Default is 50% of the entries
        Int64(floor(length(genomes.entries) * 0.5))
    else
        if (threshold_n < 1) || (threshold_n > length(genomes.entries))
            throw(
                ArgumentError(
                    "The threshold for the number of entries, `threshold_n` is out of bounds. Expected range: from 1 to " *
                    string(length(genomes.entries)) *
                    " while the supplied value is " *
                    string(threshold_n) *
                    ".",
                ),
            )
        end
        threshold_n
    end
    threshold_p = if isnothing(threshold_p)
        # Default is 50% of the traits or 1
        length(genomes.loci_alleles) == 1 ? 1 : Int64(floor(length(genomes.loci_alleles) * 0.5))
    else
        if (threshold_p < 1) || (threshold_p > length(genomes.loci_alleles))
            throw(
                ArgumentError(
                    "The threshold for the number of traits, `threshold_p` is out of bounds. Expected range: from 1 to " *
                    string(length(genomes.loci_alleles)) *
                    " while the supplied value is " *
                    string(threshold_p) *
                    ".",
                ),
            )
        end
        threshold_p
    end
    # Extract PCs

    # Impute???

    # df = tabularise(genomes)
    # traits = sort(unique(genomes.loci_alleles))
    # populations = sort(unique(genomes.populations))
    # for population in populations
    #     if sum(df.populations .== population) < 2
    #         throw(ArgumentError("Population $population has less than 2 entries."))
    #     end
    # end
    # # Output
    # (df, traits, threshold_n, threshold_p)

end
