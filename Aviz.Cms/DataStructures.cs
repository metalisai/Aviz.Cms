using System.Runtime.CompilerServices;

namespace Aviz.Cms;

using Debug = System.Diagnostics.Debug;
using V3 = System.Numerics.Vector3;

public struct EdgeCoord {
    // all this packing seems crazy, but the whole algorithms is ~10% faster
    //   with this packing
    public ulong xyz;

    public int x {
        get => (int)(xyz & 0xFFFF);
        set {
            Debug.Assert(value <= 0xFFFF);
            xyz = (xyz&0xFFFFFFFFFFFF0000UL) | ((ulong)(uint)value);
        }
    }

    public int y {
        get => (int)((xyz>>16) & 0xFFFF);
        set {
            Debug.Assert(value <= 0xFFFF);
            xyz = (xyz&0xFFFFFFFF0000FFFF) | (((ulong)(uint)value) << 16);
        }
    }

    public int z {
        get => (int)((xyz>>32) & 0xFFFF);
        set {
            Debug.Assert(value <= 0xFFFF);
            xyz = (xyz&0xFFFF0000FFFFFFFF) | (((ulong)(uint)value) << 32);
        }
    }

    // this puts hermite grid size limit at 16k
    // 16K^3 grid would take 4TB of memory, so probably ok?
    public int count {
        get => (int)((xyz>>48) & 0x3FFF);
        set {
            Debug.Assert(value <= 0x3FFF);
            xyz = (xyz&0xC000FFFFFFFFFFFF) | (((ulong)(uint)value) << 48);
        }
    }

    public int dir {
        get => (int)((xyz>>62) & 0x3);
        set {
            xyz = (xyz&0x3FFFFFFFFFFFFFFF) | (((ulong)(uint)value) << 62);
        }
    }

    public EdgeCoord() {
    }

    public EdgeCoord(I3 loc, int count, int dir) {
        this.x = loc.x;
        this.y = loc.y;
        this.z = loc.z;
        this.count = count;
        this.dir = (byte)dir;
    }

    public (EdgeCoord, EdgeCoord) Subdivide() {
        Debug.Assert(count > 0 && (count&1) == 0);
        var hc = count>>1;
        var first = new EdgeCoord(){ x = x, y = y, z = z, count = hc, dir=dir };
        var second = dir switch {
            0 => new EdgeCoord(){ x = x + hc, y = y, z = z, count = hc, dir=0 },
            1 => new EdgeCoord(){ x = x, y = y + hc, z = z, count = hc, dir=1 },
            2 or _ => new EdgeCoord(){ x = x, y = y, z = z + hc, count = hc, dir=2 },
        };
        return (first, second);
    }

    public static bool operator==(EdgeCoord a, EdgeCoord b) => a.xyz == b.xyz;
    public static bool operator!=(EdgeCoord a, EdgeCoord b) => a.xyz != b.xyz; 
    public override bool Equals(object? obj) => obj is EdgeCoord other && Equals(other);
    public bool Equals(EdgeCoord other) => xyz == other.xyz;
    public override int GetHashCode() => HashCode.Combine(xyz);
}

public abstract class Vertex {
    public abstract override bool Equals(object? obj);
    public abstract override int GetHashCode();
}

public class EdgeVertex : Vertex {
    public EdgeCoord edge { get; set; }
    public EdgeVertex(EdgeCoord edge) {
        this.edge = edge;
    }

    public override bool Equals(object? obj) => obj is EdgeVertex ev ? this.edge == ev.edge : false;
    public override int GetHashCode() => edge.GetHashCode();
}

public class NewVertex : Vertex {
    public V3 position;
    public I3 cellMin;
    public int cellSize;

    public NewVertex(V3 position, I3 cellMin, int cellSize)
    {
        this.position = position;
        this.cellMin = cellMin;
        this.cellSize = cellSize;
    }

    public override bool Equals(object? obj) => ReferenceEquals(this, obj);
    public override int GetHashCode() => System.Runtime.CompilerServices.RuntimeHelpers.GetHashCode(this);
}

public record struct Segment(Vertex v1, Vertex v2, int nDir);

[InlineArray(4)]
public struct FaceEdges {
    public EdgeCoord element0;
}

public struct Face {
    public int firstChild;
    public bool IsLeaf => firstChild == -1;
    public Face() {
        firstChild = -1;
    }
}

[InlineArray(6)]
public struct CellFaces {
    public int element0;
}

public struct Cell {
    // 8 byte reference or 64 byte array
    public int firstChild;
    public bool IsLeaf => firstChild == -1;

    // subdivision constructor
    public Cell() {
        firstChild = -1;
    }      
}
