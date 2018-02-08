using Distributions
using DataFrames

"""
    generateBootMatrix(data, factorLevels, numBootSamples)

This function takes two vector inputs - first, the data to bootstrap from and
second, a vector of factor identifiers whic represent the factor level combinations
that generated the data.

The output is a bootstrapped matrix with first column filled with the original
data and the subsequent columns filled with stratified samples of the data
"""
function generateBootMatrix(data::Vector{T},
                            factorIdentifier::Vector{T},
                            numBootSamples::Int) where T <: Real

    # Allocate the matrix for bootstrapped data
    bootMatrix = Matrix{Float64}(length(data), numBootSamples + 1)

    # Fill the first column with the sample
    bootMatrix[:, 1] = data

    for s in 2:(numBootSamples+1)
        for label in unique(factorIdentifier)
            subsetIndices = find(x -> (x == label), factorIdentifier)
            dataSubset = data[subsetIndices]
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
