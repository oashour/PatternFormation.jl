function init_cond(type, N; dx=1/143, M = 5, α=500.0, β=1000.0, c01 = 1, BC=:Periodic, r₀ = 0.0)
    # Use odd sized pattern
    @assert M ÷ 2 + 1 == ceil(M/2)
    # Done
    if type==:NoisePatches
        u01 = c01*ones(Float64, N,N)
        n_rect = rand(20:50)
        x = collect(-N/2*dx:dx:(N/2-1)*dx) # Global coordinate
        dy = dx*sqrt(3)/2
        y = collect(-N/2*dy:dy:(N/2-1)*dy) # Global Coordinate
        for i = 1:n_rect 
            l = rand(1:N ÷ 5)
            w = rand(1:N ÷ 5)
            sx = rand(1:4*N ÷ 5)
            sy = rand(1:4*N ÷ 5)
            val = rand()
            u01[sx:sx+l,sy:sy+l] .= val
        end
    elseif type==:Square
        @assert M ≥ 3
        x = collect(-N/2*dx:dx:(N/2-1)*dx) # Global coordinate
        y = collect(-N/2*dx:dx:(N/2-1)*dx) # Global Coordinate
        len = x[end]/(M-1)*2 # The factor here is to make the whole MxM grid fit into the lattice
        basis = len*hcat([0,1],[1,0]) |> SMatrix{2,2}
        cell = HomogeneousCell([[0.0,0.0]])
        lattice = RegularLattice((M,M), basis, cell)
        u01 = zeros(Float64, N,N)
        for i in 1:M
            for j in 1:M
                ii, jj= (lattice[i,j,1] .- len*(M+1)/2)
                u01 .+= exp.(-α*((x .- ii).^2 .+ (y' .- jj).^2))
            end
        end
    elseif type==:ErMnO3
        @assert M ≥ 5
        # Grid Stuff
        x = collect(-N/2*dx:dx:(N/2-1)*dx) # Global coordinate
        dy = dx*sqrt(3)/2 # factor of 6/5??
        y = collect(-N/2*dy:dy:(N/2-1)*dy) # Global Coordinate
        #y = [y; y[end]*(1+dx)] # Hack Fix, idk why y sometimes has N-1 elements
        u01 = zeros(Float64, N,N)
        len = x[end]/(M-3)*2
        # Declare Basis
        fbasis = len*hcat([0.5, 0.5*sqrt(3)],
                      [0.5, -0.5*sqrt(3)],
                     ) |> SMatrix{2,2}
        # Atom positions
        cell_vectors_raw1 = [[0.0, 0.0], [1/3, 2/3], [2/3, 1/3]]
        cell_vectors_raw2 = [[0, 1/3], [0, 2/3], [1/3, 1/3], [1/3, 0], [2/3, 0], [2/3, 2/3]]

        fcell = InhomogeneousCell([fbasis*vec for vec in cell_vectors_raw1], [fbasis*vec for vec in cell_vectors_raw2]; label = :ermno3)
        lattice = RegularLattice((M,M), fbasis, fcell; label = :hexagonal)

        u01 = zeros(Float64, N,N)
        for i in 1:M
            for j in 1:M
                for ic in 1:length(cell_vectors_raw1)
                    ii, jj= lattice[i,j,ic,1]
                    ii -= len*(M ÷ 2 + 1)
                    if BC != :Periodic
                      if isapprox(ii,x[end],atol=dx) || isapprox(ii,x[1],atol=dx) || isapprox(jj,y[1],atol=dx) || isapprox(jj,y[end],atol=dx)
                        #println("on boundary")
                        continue
                      else
                        x₀ = r₀*x[end]*rand()
                        y₀ = r₀*x[end]*rand()
                        u01 .+= exp.(-α*((x .- ii .- x₀).^2 .+ (y' .- jj .- y₀).^2))
                      end
                    else
                        u01 .+= exp.(-α*((x .- ii).^2 .+ (y' .- jj).^2))
                    end
                end
                for ic in 1:length(cell_vectors_raw2)
                    ii, jj= lattice[i,j,ic,2]
                    ii -= len*(M ÷ 2 + 1)
                    if BC != :Periodic
                      if isapprox(ii,x[end],atol=dx) || isapprox(ii,x[1],atol=dx) || isapprox(jj,y[1],atol=dx) || isapprox(jj,y[end],atol=dx)
                        #println("on boundary")
                        continue
                      else
                        x₀ = r₀*x[end]*rand()
                        y₀ = r₀*x[end]*rand()
                        u01 .+= exp.(-β*((x .- ii .- x₀).^2 .+ (y' .- jj .- y₀).^2))
                      end
                    else
                        u01 .+= exp.(-β*((x .- ii).^2 .+ (y' .- jj).^2))
                    end
                end
            end
        end
    elseif type==:ErMnO3_ErOnly
        @assert M ≥ 5
        # Grid Stuff
        x = collect(-N/2*dx:dx:(N/2-1)*dx) # Global coordinate
        dy = dx*sqrt(3)/2 # factor of 6/5??
        y = collect(-N/2*dy:dy:(N/2-1)*dy) # Global Coordinate
        #y = [y; y[end]*(1+dx)] # Hack Fix, idk why y sometimes has N-1 elements
        u01 = zeros(Float64, N,N)
        len = x[end]/(M-3)*2
        # Declare Basis
        fbasis = len*hcat([0.5, 0.5*sqrt(3)],
                      [0.5, -0.5*sqrt(3)],
                     ) |> SMatrix{2,2}
        # Atom positions
        cell_vectors_raw1 = [[0.0, 0.0], [1/3, 2/3], [2/3, 1/3]]
        cell_vectors_raw2 = [[0, 1/3], [0, 2/3], [1/3, 1/3], [1/3, 0], [2/3, 0], [2/3, 2/3]]

        fcell = InhomogeneousCell([fbasis*vec for vec in cell_vectors_raw1], [fbasis*vec for vec in cell_vectors_raw2]; label = :ermno3)
        lattice = RegularLattice((M,M), fbasis, fcell; label = :hexagonal)

        u01 = zeros(Float64, N,N)
        for i in 1:M
            for j in 1:M
                for ic in 1:length(cell_vectors_raw1)
                    ii, jj= lattice[i,j,ic,1]
                    ii -= len*(M ÷ 2 + 1)
                    if BC != :Periodic
                      if isapprox(ii,x[end],atol=dx) || isapprox(ii,x[1],atol=dx) || isapprox(jj,y[1],atol=dx) || isapprox(jj,y[end],atol=dx)
                        #println("On Boundary")
                        continue
                      else
                        x₀ = r₀*x[end]*rand()
                        y₀ = r₀*x[end]*rand()
                        u01 .+= exp.(-α*((x .- ii .- x₀).^2 .+ (y' .- jj .- y₀).^2))
                      end
                      else
                        u01 .+= exp.(-α*((x .- ii).^2 .+ (y' .- jj).^2))
                    end
                end
            end
        end
    elseif type==:ErMnO3_MnOnly
        @assert M ≥ 5
        # Grid Stuff
        x = collect(-N/2*dx:dx:(N/2-1)*dx) # Global coordinate
        dy = dx*sqrt(3)/2 # factor of 6/5??
        y = collect(-N/2*dy:dy:(N/2-1)*dy) # Global Coordinate
        #y = [y; y[end]*(1+dx)] # Hack Fix, idk why y sometimes has N-1 elements
        u01 = zeros(Float64, N,N)
        len = x[end]/(M-3)*2
        # Declare Basis
        fbasis = len*hcat([0.5, 0.5*sqrt(3)],
                      [0.5, -0.5*sqrt(3)],
                     ) |> SMatrix{2,2}
        # Atom positions
        cell_vectors_raw1 = [[0.0, 0.0], [1/3, 2/3], [2/3, 1/3]]
        cell_vectors_raw2 = [[0, 1/3], [0, 2/3], [1/3, 1/3], [1/3, 0], [2/3, 0], [2/3, 2/3]]

        fcell = InhomogeneousCell([fbasis*vec for vec in cell_vectors_raw1], [fbasis*vec for vec in cell_vectors_raw2]; label = :ermno3)
        lattice = RegularLattice((M,M), fbasis, fcell; label = :hexagonal)

        u01 = zeros(Float64, N,N)
        for i in 1:M
            for j in 1:M
                for ic in 1:length(cell_vectors_raw2)
                    ii, jj= lattice[i,j,ic,2]
                    ii -= len*(M ÷ 2 + 1)
                    if BC != :Periodic
                      if isapprox(ii,x[end],atol=dx) || isapprox(ii,x[1],atol=dx) || isapprox(jj,y[1],atol=dx) || isapprox(jj,y[end],atol=dx)
                        continue
                      else
                        x₀ = r₀*x[end]*rand()
                        y₀ = r₀*x[end]*rand()
                        u01 .+= exp.(-β*((x .- ii .- x₀).^2 .+ (y' .- jj .- y₀).^2))
                      end
                    else
                        u01 .+= exp.(-β*((x .- ii).^2 .+ (y' .- jj).^2))
                    end
                end
            end
        end
    elseif type==:HoneyComb
        @assert M ≥ 5
        # Grid Stuff
        x = collect(-N/2*dx:dx:(N/2-1)*dx) # Global coordinate
        dy = dx*sqrt(3)/2*6/5 # factor of 6/5
        y = collect(-N/2*dy:dy:(N/2-1)*dy) # Global Coordinate
        #y = [y; y[end]*(1+dx)] # Hack Fix, idk why y sometimes has N-1 elements
        u01 = zeros(Float64, N,N)
        len = x[end]/(M-3)*2
        # Declare Basis
        fbasis = len*hcat([0.5, 0.5*sqrt(3)],
                      [0.5, -0.5*sqrt(3)],
                     ) |> SMatrix{2,2}
        # Atom positions
        cell_vectors_raw = [[2/3, 1/3], [1/3, 2/3]]
        pos = [fbasis*vec for vec in cell_vectors_raw]
        fcell = HomogeneousCell(pos)
        lattice = RegularLattice((M,M), fbasis, fcell; label = :hexagonal)

        u01 = zeros(Float64, N,N)
        for i in 1:M
            for j in 1:M
                for ic in 1:length(cell_vectors_raw)
                    ii, jj= lattice[i,j,ic]
                    ii -= len*(M ÷ 2 + 1)
                    u01 .+= exp.(-α*((x .- ii).^2 .+ (y' .- jj).^2))
                end
            end
        end
    elseif type==:Hexagonal
        @assert M ≥ 5
        # Grid Stuff
        x = collect(-N/2*dx:dx:(N/2-1)*dx) # Global coordinate
        dy = dx*sqrt(3)/2
        y = collect(-N/2*dy:dy:(N/2-1)*dy) # Global Coordinate
        #y = [y; y[end]*(1+dx)] # Hack Fix, idk why y sometimes has N-1 elements
        u01 = zeros(Float64, N,N)
        len = x[end]/(M-3)*2
        # Declare Basis
        fbasis = len*hcat([0.5, 0.5*sqrt(3)],
                      [0.5, -0.5*sqrt(3)],
                     ) |> SMatrix{2,2}
        # Atom positions
        cell_vectors_raw = [[0,0]]
        pos = [fbasis*vec for vec in cell_vectors_raw]
        fcell = HomogeneousCell(pos)
        lattice = RegularLattice((M,M), fbasis, fcell; label = :hexagonal)

        for i in 1:M
            for j in 1:M
                for ic in 1:length(cell_vectors_raw)
                    ii, jj= lattice[i,j, ic]
                    ii -= len*(M ÷ 2 + 1)
                    u01 .+= exp.(-α*((x .- ii).^2 .+ (y' .- jj).^2))
                end
            end
        end
    end
    return u01, x, y, dx, dy
end
