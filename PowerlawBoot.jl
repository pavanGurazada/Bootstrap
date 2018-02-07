using Distributions
srand(20130810)

# Generate a random sample from a Pareto distribution
numSamples = 10^4
data = rand(Pareto(), numSamples)

function generateBootMatrix(data::Vector{T}, numBootSamples::Int) where T <: Real
    # Allocate the boot matrix
    bootMatrix = Matrix{Float64}(length(data), numBootSamples)

    # Fill the first column with the sample
    bootMatrix[:, 1] = data

    # Fill the rest of the columns with the bootstrapped data

    for s in 2:numBootSamples
        bootMatrix[:, s] = sample(data, numSamples, replace = true)
    end

    return bootMatrix
end
