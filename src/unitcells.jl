# TODO: ほしい格子に対して unit cell データを追加していく
"""
    get_unit_cell(::Type{T}) where T <: AbstractTopology
各トポロジーに対応する単位胞データを返す関数。
"""
function get_unit_cell(::Type{T}) where T<:AbstractTopology
    error("UnitCell not defined for $T")
end

struct Square <: AbstractTopology{2} end
function get_unit_cell(::Type{Square})
    a1 = [1.0, 0.0]
    a2 = [0.0, 1.0]
    conns = [
        Connection(1, 1, 1, 0, 1),
        Connection(1, 1, 0, 1, 2)
    ]
    return UnitCell{2,Float64}([a1, a2], [[0.0, 0.0]], conns)
end

struct Triangular <: AbstractTopology{2} end
function get_unit_cell(::Type{Triangular})
    a1 = [1.0, 0.0]
    a2 = [0.5, sqrt(3) / 2]
    conns = [
        Connection(1, 1, 1, 0, 1), # 右 (a1)
        Connection(1, 1, 0, 1, 2), # 上 (a2)
        Connection(1, 1, -1, 1, 3) # 左上 (a2-a1)
    ]
    return UnitCell{2,Float64}([a1, a2], [[0.0, 0.0]], conns)
end

struct Honeycomb <: AbstractTopology{2} end
function get_unit_cell(::Type{Honeycomb})
    a1 = [sqrt(3), 0.0]
    a2 = [sqrt(3) / 2, 1.5]
    d_A = [0.0, 0.0]
    d_B = [0.0, 1.0]
    conns = [
        Connection(1, 2, 0, 0, 1),  # A -> B (同セル内: 上)
        Connection(1, 2, 1, -1, 2), # A -> B (左上のセル)
        Connection(1, 2, 0, -1, 3)   # A -> B (上のセル)
    ]
    return UnitCell{2,Float64}([a1, a2], [d_A, d_B], conns)
end

struct Kagome <: AbstractTopology{2} end
function get_unit_cell(::Type{Kagome})
    a1 = [1.0, 0.0]
    a2 = [0.5, sqrt(3) / 2]
    d_A = [0.0, 0.0]
    d_B = 0.5 * a1
    d_C = 0.5 * a2
    conns = [
        # 単位胞内の三角形
        Connection(1, 2, 0, 0, 1), # A-B
        Connection(1, 3, 0, 0, 1), # A-C
        Connection(2, 3, 0, 0, 1), # B-C
        # 隣の単位胞への接続
        Connection(2, 1, 1, 0, 1), # B -> Next A (右)
        Connection(3, 1, 0, 1, 1), # C -> Next A (上)
        Connection(2, 3, 1, -1, 1) # B -> Next C (右下: Kagome特有のクロス)
    ]
    return UnitCell{2,Float64}([a1, a2], [d_A, d_B, d_C], conns)
end
struct Lieb <: AbstractTopology{2} end
function get_unit_cell(::Type{Lieb})
    # 正方格子の辺中心にサイトを追加した格子 (CuO2面など)
    a1 = [2.0, 0.0]
    a2 = [0.0, 2.0]

    # A:Corner, B:Right, C:Up
    d_A = [0.0, 0.0]
    d_B = [1.0, 0.0]
    d_C = [0.0, 1.0]

    conns = [
        Connection(1, 2, 0, 0, 1), # A -> B (Intra)
        Connection(1, 3, 0, 0, 2), # A -> C (Intra)
        Connection(2, 1, 1, 0, 1), # B -> Next A (Right)
        Connection(3, 1, 0, 1, 2)  # C -> Next A (Up)
    ]
    return UnitCell{2,Float64}([a1, a2], [d_A, d_B, d_C], conns)
end
const AVAILABLE_LATTICES = (Square, Triangular, Honeycomb, Kagome, Lieb)
export AVAILABLE_LATTICES
export Square, Triangular, Honeycomb, Kagome, Lieb