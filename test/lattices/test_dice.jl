@testset "Dice Lattice Tests" begin
    @testset "Unit Cell Definition" begin
        uc = get_unit_cell(Dice)
        # Basis Vectors are separated by 60 degrees and unit length (triangular)
        a1, a2 = uc.basis
        @test dot(a1, a2) ≈ cos(π/3) atol=1e-10
        @test norm(a1) ≈ sqrt(1) atol=1e-10
        @test norm(a2) ≈ sqrt(1) atol=1e-10
        # 3 sublattices: Hub(1), RimA(2), RimB(3)
        @test length(uc.sublattice_positions) == 3

        d_1 = uc.sublattice_positions[1]
        d_2 = uc.sublattice_positions[2]
        d_3 = uc.sublattice_positions[3]
        @test d_1 == [0.0, 0.0]
        @test d_2 ≈ (a1 .+ a2) ./ 3
        @test d_3 ≈ (a1 .+ a2) .* (2 / 3)
    end

    @testset "Geometry & Reciprocal Vectors" begin
        Lx, Ly = 4, 4
        lat = build_lattice(Dice, Lx, Ly)

        # a is the basis vectors, b is the reciprocal vectors
        # Verify a_i · b_j = 2π δ_ij
        a = lat.basis_vectors
        b = lat.reciprocal_vectors

        @test dot(a[1], b[1]) ≈ 2π atol=1e-10
        @test dot(a[2], b[2]) ≈ 2π atol=1e-10
        @test abs(dot(a[1], b[2])) < 1e-10
        @test abs(dot(a[2], b[1])) < 1e-10
    end

    @testset "Topology & Connectivity" begin
        Lx, Ly = 4, 4
        @testset "periodic boundary condition" begin
            lat_pbc = build_lattice(Dice, Lx, Ly; boundary=PBC())

            # total number of sites will be 3 * Lx * Ly
            @test lat_pbc.N == 3 * Lx * Ly
            
            # Dice lattice: Hub has coordination 6, Rims have coordination 3
            degrees = length.(lat_pbc.nearest_neighbors)
            sub_ids = lat_pbc.sublattice_ids
            
            # Hub (sublattice 1) should have degree 6
            @test all(degrees[sub_ids .== 1] .== 6)
            # Rim A and B (sublattices 2, 3) should have degree 3
            @test all(degrees[sub_ids .== 2] .== 3)
            @test all(degrees[sub_ids .== 3] .== 3)
            
            id_Hub = lat_pbc.site_map[1, 1]
            id_RimA = id_Hub + 1
            id_RimB = id_Hub + 2
            # Check if the neighbors contain the corresponding IDs
            @test id_RimA in lat_pbc.nearest_neighbors[id_Hub]
            @test id_RimB in lat_pbc.nearest_neighbors[id_Hub]
            @test id_Hub in lat_pbc.nearest_neighbors[id_RimA]

            # Bipartite check - Dice is bipartite (Hub vs Rim)
            @test lat_pbc.is_bipartite == true
        end
        @testset "open boundary condition" begin
            lat_obc = build_lattice(Dice, Lx, Ly; boundary=OBC())
            
            degrees = length.(lat_obc.nearest_neighbors)
            @test maximum(degrees) <= 6
            @test minimum(degrees) >= 1

            # Bulk Hub site should have coordination close to 6
            # Bulk Rim site should have coordination close to 3
            bulk_hub = lat_obc.site_map[2, 2]
            @test length(lat_obc.nearest_neighbors[bulk_hub]) <= 6
        end
    end
    @testset "Index Consistency" begin
        Lx, Ly = 3, 3
        n_sub = 3 # Dice has 3 sublattices
        lat = build_lattice(Dice, Lx, Ly)
        uc = get_unit_cell(Dice)
        
        for x in 1:Lx, y in 1:Ly
            base_idx = lat.site_map[x, y]
            
            # Base index calc check
            expected_base = ((x - 1) + (y - 1) * Lx) * n_sub + 1
            @test base_idx == expected_base
            
            # Sublattice position check
            cell_origin = (x-1)*lat.basis_vectors[1] + (y-1)*lat.basis_vectors[2]
            
            # Check Hub, RimA, RimB positions
            for s in 1:n_sub
                pos = lat.positions[base_idx + s - 1]
                expected_pos = cell_origin + uc.sublattice_positions[s]
                @test pos ≈ expected_pos
            end
        end
    end
end
