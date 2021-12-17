
<<<<<<< HEAD
function set_data_to_ones(A)
	# create copy of matrix A
	matrix = copy(A)
	# grab dimension
	dim = size(matrix)[1]
	# loop thru matrix, changing all non zero entries into 1's
	for i = 1 : dim
		for j = 1 : dim
			if matrix[i, j] != 0
				matrix[i, j] = 1
			end
		end
	end
	return matrix
end

function degree_matrix(A)
	# grab matrix dimension
	dim = size(A)[1]
	# allocate diagonal vector
	d = zeros(dim)
	# loop thru matrix
	for i = 1 : dim
		for k = 1 : dim
			# set diagonal vector 'd' to the sum of non zero values in cooresponding row of 'A'
			d[i] = d[i] + A[i, k]
		end
	end
	# turn 'd' into sparse diagonal representation
	D = sparse(Diagonal(d))
	return D
end

function create_matrices(A)
	temp = findnz(A)
	row = temp[1]
	col = temp[2]
	data = ones(nnz(A))
	Aₐ = sparse(row, col, data)
	w = rand(nnz(Aₐ))
	dim = size(Aₐ)[1]
	E = spzeros(nnz(Aₐ), dim)
		for i=1:nnz(Aₐ)
			first = row[i]
			second = col[i]
			E[i, first] = 1
			E[i, second] = 1
		end
	Eᵀ = sparse(E')
	EEᵀ = E * Eᵀ
	return Aₐ, E, Eᵀ, EEᵀ, w
end

# multiples list of intermediate coarsened 'P' matrices, returning one matrix
function interpolate_P(Pᵢ)
	k = length(Pᵢ)
	P = Pᵢ[1]
	for i = 2 : k
		P = P * Pᵢ[i]
	end
	Pᵀ = sparse(P')
	return P, Pᵀ
end

function luby(A, EEᵀ, w)
	label  = (-1).* ones(length(w))
	matching_counter = 1
	for i=1:length(w)
		edge_is_max = true
		related_edges = findnz(EEᵀ[i,:])
		for k in related_edges[1]
			if w[i] < w[k]
			 	edge_is_max = false
			end
		end
		if edge_is_max == true
			label[i] = matching_counter
			matching_counter = matching_counter + 1
		end
	return matching_counter, label
end

function create_aggragates(matching_counter, Eᵀ, label)
	aggragate_counter = matching_counter
	vertex = size(Eᵀ)[1]
	vertex_aggragate = zeros(vertex)
	for i = 1:vertex
		not_shared = true
		related_edges = findnz(Eᵀ[i,:])[1]
		for k in related_edges
			if label[k] > -1
				vertex_aggragate[i] = label[k]
				not_shared = false
			end
		end
		if not_shared == true
			vertex_aggragate[i] = aggragate_counter
			aggragate_counter = aggragate_counter + 1
		end
	end
	return vertex_aggragate
end

function create_P(vertex_aggragate)
	dim_n = maximum(vertex_aggragate)
	dim_n = convert(Int, dim_n)
	dim_m = length(vertex_aggragate)[1]
	P = spzeros(dim_m, dim_n)
	for i=1:dim_m
		position = convert(Int, vertex_aggragate[i])
		P[i, position] = 1
	end
	Pᵀ = sparse(P')
	return P, Pᵀ
end

function recurse(A, EEᵀ, Eᵀ, w, τ, A_nnz, Aᵢ, Pᵢ)
	matching_counter, label = luby(A, EEᵀ, w)
	vertex_aggragate = create_aggragates(matching_counter, Eᵀ, label)
	P, Pᵀ = create_P(vertex_aggragate)
	push!(Aᵢ, A)
	push!(Pᵢ, P)
	PᵀAP = Pᵀ * A * P
	if nnz(PᵀAP) <= A_nnz/τ
		return PᵀAP, Aᵢ, Pᵢ
	else
		A, E, Eᵀ, EEᵀ, w = create_matrices(PᵀAP)
		recurse(A, EEᵀ, Eᵀ, w, τ, A_nnz, Aᵢ, Pᵢ)
	end
end

function forward_sidel(A)
    return lowerTri=tril(A)
end

function backward_sidel(A)
    return upperTri=triu(A)
end

function two_level(r, A, PᵀAP, P, Pᵀ)
	M = forward_sidel(A)
	Mᵀ = backward_sidel(A)
	y = M \ r
	rᶜ = Pᵀ * (r - A * y)
	yᶜ = PᵀAP \ rᶜ
	y = y + (P * yᶜ)
	z = Mᵀ \ (r - A * y)
	y = y + z
	r̄ = y
	return r̄
end

function l1_smoother(A)
	dim = size(A)[1]
	l =  zeros(dim)
	for i = 1 : dim
		for k = 1 : dim
			l[i] = l[i] + abs(A[i, k])
		end
	end
	L1 = sparse(Diagonal(l))
	return L1
end

function B_matrix(A, Aᶜᵢ, Pᵢ, Pᵀᵢ, r, k)
	if k == 1
		l1 = l1_smoother(A)
		r̄ = l1 \ r
		return r̄
	else
		r̄ = two_level(r, A, Aᶜᵢ[k-1], Pᵢ[k-1], Pᵀᵢ[k-1])
		return r̄
	end
end

function composite_solver(A, b, x, iterations, Aᶜᵢ, Pᵢ, Pᵀᵢ, τ)
	k = length(Aᶜᵢ) + 1
	y = zeros(size(x)[1])
	r = b - A * x
	for i = 1 : k
		y = B_matrix(A, Aᶜᵢ, Pᵢ, Pᵀᵢ, r, i)
		x = x + y
		r = r - A * y
	end
	for i = k : -1 : 1
		y = B_matrix(A, Aᶜᵢ, Pᵢ, Pᵀᵢ, r, i)
		x = x + y
		r = r - A * y
	end
	if iterations > 1
		iterations -= 1
		composite_solver(A, b, x, iterations, Aᶜᵢ, Pᵢ, Pᵀᵢ, τ)
	end
	return x
end

function wrapper(A, b, x, Aᶜᵢ, Pᵢ, Pᵀᵢ, Wᵢ, τ, k, iterations)
	xᵢ = composite_solver(A, b, x, iterations, Aᶜᵢ, Pᵢ, Pᵀᵢ, τ)
	wᵢ = (1 / norm(xᵢ)) .* xᵢ
	push!(Wᵢ, wᵢ)
	if k == 1
		return Wᵢ
	end
	Aₐ, E, Eᵀ, EEᵀ, rand = create_matrices(A)
	A_nnz = nnz(A)
	Aₖ = [ ]
	Pₖ = [ ]
	Āᶜ, Aₖ, Pₖ = recurse(A, EEᵀ, Eᵀ, rand, τ, A_nnz, Aₖ, Pₖ)
	P̄, P̄ᵀ = interpolate_P(Pₖ)
	W = sparse(Diagonal(wᵢ))
	P = W * P̄
	push!(Pᵢ, P)
	Pᵀ = sparse(P')
	push!(Pᵀᵢ, Pᵀ)
	Aᶜ = Pᵀ * A * P
	push!(Aᶜᵢ, Aᶜ)
	k -= 1
	wrapper(A, b, x, Aᶜᵢ, Pᵢ, Pᵀᵢ, Wᵢ, τ, k, iterations)
end

end
end
=======
>>>>>>> 1a5ff369269e5d607a625b2400c86c0487fb703d
