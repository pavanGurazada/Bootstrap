using Distributions
using DataFrames

"""
    generateBootMatrix(y, [X, numBootSamples])

This function takes two vector inputs - first, the data to bootstrap from and
second, a vector of factor identifiers which represent the factor level combinations
that generated the data.

The output is a bootstrapped matrix with first column filled with the original
data and the subsequent columns filled with stratified samples of the data.

The function is set up in this way for minimal dependence on other packages like
DataFrames that are still in a very primitive stage of development. It is easy
to get to the required inputs of this function from a DataFrame-ish sort of an
object by extracting the vector of outcome variable (y) and a vector of factor
level combination identifiers (X).
"""
function generateBootMatrix(y::Vector{T},
                            X::Vector{T},
                            numBootSamples::Int) where T <: Real

    # Allocate the matrix for bootstrapped data
    bootMatrix = Matrix{Float64}(length(y), numBootSamples + 1)

    # Fill the first column with the sample
    bootMatrix[:, 1] = y

    for s in 2:(numBootSamples+1)
        for label in unique(X)
            subsetIndices = find(x -> (x == label), X)
            dataSubset = y[subsetIndices]
            bootMatrix[subsetIndices, s] = sample(dataSubset, length(dataSubset), replace = true)
        end
    end

    return bootMatrix
end

# Minimal working examples

# Generate a random sample
numSamplesPerFactor = 100

data = append!(rand(Pareto(), numSamplesPerFactor),
               rand(Uniform(), numSamplesPerFactor))
factorIdentifier = append!(ones(numSamplesPerFactor),
                           zeros(numSamplesPerFactor))

numBootSamples = 10^4

# Perform the bootstrap
@time bootData = generateBootMatrix(data, factorIdentifier, numBootSamples)
