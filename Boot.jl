using Distributions
srand(20130810)

"""
    generateBootMatrix(data, n)

This is a general purpose function that builds a matrix with bootstrapped samples
from the data. The first column of the output is the original data. Subsequent columns
contain the bootstrapped samples.
"""
function generateBootMatrix(data::Vector{T}, numBootSamples::Int) where T <: Real
    # Allocate the boot matrix
    bootMatrix = Matrix{Float64}(length(data), numBootSamples+1)

    # Fill the first column with the sample
    bootMatrix[:, 1] = data

    # Fill the rest of the columns with the bootstrapped data
    for s in 2:(numBootSamples+1)
        bootMatrix[:, s] = sample(data, numSamples, replace = true)
    end

    return bootMatrix
end

"""
    summarize(m[, f])

Summarize the columns of a given boot matrix with the provided function.

...
# Arguments
- `m::Matrix`: the bootstrapped matrix
- `f::Function`: the function to apply on the columns
...

First, allocate the output and apply `f` to each of the
columns. Return the output in a "tall form" with the original sample summary as the
first row and the subsequent rows containing the summarized data from the other
columns of the matrix.

Supported summarizing functions are those that take a vector and return a
real (e.g., mean, median)
"""
function summarize(bootMatrix::Array{Float64, 2}, f::Function)
    ncols = size(bootMatrix)[2]

    result = Matrix{Float64}(ncols, 2)
    result[:, 1] = 1:ncols

    result[:, 2] = mapslices(f, bootMatrix[:, :], 1)

    return result
end

# Minimal working examples

# Generate a random sample from a Pareto distribution
numSamples = 100
data = rand(Pareto(), numSamples)
numBootSamples = 10^4

@time bootData = generateBootMatrix(data, numBootSamples)
@time bootMeans = summarize(bootData, mean)
@time bootMedians = summarize(bootData, median)
