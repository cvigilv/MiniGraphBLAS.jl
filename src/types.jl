# Matrix
struct GBMatrix{T} <: AbstractMatrix{T}
    data::Matrix{T}
end

GBMatrix(data::Matrix{T}) where T = GBMatrix{T}(data)
GBMatrix(rows::Int, cols::Int) = GBMatrix(zeros(rows, cols))

Base.size(A::GBMatrix) = size(A.data)
Base.getindex(A::GBMatrix, i::Int, j::Int) = A.data[i, j]
Base.setindex!(A::GBMatrix, v, i::Int, j::Int) = (A.data[i, j] = v)


# Vector
struct GBVector{T} <: AbstractVector{T}
    data::Vector{T}
end

GBVector(data::Vector{T}) where T = GBVector{T}(data)
GBVector(size::Int) = GBVector(zeros(size))

Base.size(A::GBVector) = size(A.data)
Base.getindex(A::GBVector, i::Int) = A.data[i]
Base.setindex!(A::GBVector, v, i::Int) = (A.data[i] = v)


# Linear algebra primitives
"""
    Base.:*(A::GBMatrix{T}, B::GBMatrix{T}, oplus::Function=+, otimes::Function=*)

Matrix times matrix. Process connecting outgoing edges.

# Arguments
- `A::GBMatrix{T}`: GraphBLAS matrix
- `B::GBMatrix{T}`: GraphBLAS matrix
- `oplus::Function`: Addition operator
- `otimes::Function`: Multiplication operator
"""
function Base.:*(
	A::GBMatrix{T},
	B::GBMatrix{T},
	oplus::Function=+,
	otimes::Function=*,
) where T
	# Get matrices dimensions
    rA, cA = size(A)
    rB, cB = size(B)

	# Check if multiplication is valid
    if cA != rB
        error("Matrix dimensions mismatch for multiplication")
    end

	# Initialize results matrix
	C = zeros(T, rA, cB)

	# Do matrix multiplication
    for i in 1:rA
        for j in 1:cB
            for k in 1:cA
				C[i,j] = oplus(C[i,j], otimes(A[i,k], B[k,j]))
            end
		end
    end

	return GBMatrix(C)
end

"""
    Base.:*(A::GBMatrix{T}, v::GBVector{T}, oplus::Function=+, otimes::Function=*)

Matrix times vector. Process incoming edges.

# Arguments
- `A::GBMatrix{T}`: GraphBLAS matrix
- `v::GBVector{T}`: GraphBLAS vector
- `oplus::Function`: Addition operation
- `otimes::Function`: Multiplication operation
"""
function Base.:*(
	A::GBMatrix{T},
	v::GBVector{T},
	oplus::Function=+,
	otimes::Function=*,
) where T
    m, n = size(A)
    if n != length(v)
        error("Matrix columns must match vector length")
    end

    M = zeros(T, m)

    for i in 1:m
        for j in 1:n
			M[i] = oplus(M[i], otimes(A[i,j],v[j]))
        end
    end

	return GBVector(M)
end


"""
    Base.:*(v::GBVector{T}, A::GBMatrix{T}, oplus::Function=+, otimes::Function=*)

Vector times matrix. Process outgoing edges.

# Arguments
- `v::GBVector{T}`: GraphBLAS vector
- `A::GBMatrix{T}`: GraphBLAS matrix
- `oplus::Function`: Adittion operation
- `otimes::Function`: Multiplication operation
"""
function Base.:*(
	v::GBVector{T},
	A::GBMatrix{T},
	oplus::Function=+,
	otimes::Function=*,
) where T
    n = length(v)
    m, p = size(A)

    if n != m
        error("Vector length must match number of matrix rows")
    end

    M = zeros(T, p)

    for j in 1:p
        for i in 1:n
			M[j] = oplus(M[j], otimes(v[i], A[i,j]))
        end
    end

	return GBVector(M)
end
