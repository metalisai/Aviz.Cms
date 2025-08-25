using System.Runtime.CompilerServices;

namespace Aviz.Cms;

using V3 = System.Numerics.Vector3;
using V2 = System.Numerics.Vector2;

public static class Lookups {
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static I3 CornerOffsetI(int corner, int size) {
        return new I3((corner&1)!=0?size:0, (corner&2)!=0?size:0, (corner&4)!=0?size:0);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (I3 s, I3 e) EdgeOffset(int edge, int size) {
        int axis = edge >> 2;
        int idx = edge & 0x3;
        int start = axis switch {
            0 => ((idx<<1)&2) | ((idx<<1)&4), // 0, b0, b1
            1 => ((idx>>1)&1) | ((idx<<2)&4), // b1, 0, b0
            2 or _ => idx, // b0, b1, 0
        };
        int end = start | (0x1 << axis);
        return (CornerOffsetI(start, size), CornerOffsetI(end, size));
    }

    public static readonly byte[,] faceEdgeToCubeEdge = new byte[,] {
        { 4,  8, 5, 10 },
        { 6,  9, 7, 11 },
        { 8,  0, 9,  2 },
        {10,  1,11,  3 },
        { 0,  4, 1,  6 },
        { 2,  5, 3,  7 },
    };

    public static readonly byte[,] faceCubeCorners = new byte[,] {
        { 0, 2, 4, 6 },
        { 1, 3, 5, 7 },
        { 0, 4, 1, 5 },
        { 2, 6, 3, 7 },
        { 0, 1, 2, 3 },
        { 4, 5, 6, 7 },
    };

    public static readonly int [,] quadSegments = new int[18, 4] {
        { -1, -1, -1, -1},
        {  1,  0, -1, -1},
        {  0,  3, -1, -1},
        {  1,  3, -1, -1},
        {  2,  1, -1, -1},
        {  2,  0, -1, -1},
        {  0,  3,  2,  1}, // 6a
        {  2,  3, -1, -1},
        {  3,  2, -1, -1},
        {  1,  0,  3,  2}, // 9a
        {  0,  2, -1, -1},
        {  1,  2, -1, -1},
        {  3,  1, -1, -1},
        {  3,  0, -1, -1},
        {  0,  1, -1, -1},
        { -1, -1, -1, -1},
        {  0,  1,  2,  3}, // 6b
        {  1,  2,  3,  0}, // 9b
    };

    public static readonly byte[,] childOuterFaceMap = new byte[,] {
        {0b000, 0}, {0b010, 1}, {0b100, 2}, {0b110, 3},
        {0b001, 0}, {0b011, 1}, {0b101, 2}, {0b111, 3},
        {0b000, 0}, {0b001, 2}, {0b100, 1}, {0b101, 3},
        {0b010, 0}, {0b011, 2}, {0b110, 1}, {0b111, 3},
        {0b000, 0}, {0b001, 1}, {0b010, 2}, {0b011, 3},
        {0b100, 0}, {0b101, 1}, {0b110, 2}, {0b111, 3}
    };

    public static readonly byte[,] childInnerFaceMap = new byte[,] {
        {0, 1}, {1, 0}, {2, 1}, {3, 0},
        {4, 1}, {5, 0}, {6, 1}, {7, 0},
        {0, 3}, {2, 2}, {4, 3}, {6, 2},
        {1, 3}, {3, 2}, {5, 3}, {7, 2},
        {0, 5}, {4, 4}, {1, 5}, {5, 4},
        {2, 5}, {6, 4}, {3, 5}, {7, 4},
    };

    public static readonly byte[,] faceInnerEdgeMap = new byte[,] {
        {2, 2}, {3, 4}, {2, 3}, {1, 4},
        {0, 2}, {3, 5}, {0, 3}, {1, 5},
        {2, 4}, {3, 0}, {2, 5}, {1, 0},
        {0, 4}, {3, 1}, {0, 5}, {1, 1},
        {2, 0}, {3, 2}, {2, 1}, {1, 2},
        {0, 0}, {3, 3}, {0, 1}, {1, 3},
    };

    public static readonly byte[,] faceOuterEdgeMap = new byte[,] {
        {0, 4, 2}, {1, 2, 0}, {0, 4, 3}, {3, 3, 0},
        {1, 2, 1}, {2, 5, 2}, {2, 5, 3}, {3, 3, 1},
        {0, 0, 2}, {1, 4, 0}, {0, 0, 3}, {3, 5, 0},
        {1, 4, 1}, {2, 1, 2}, {2, 1, 3}, {3, 5, 1},
        {0, 2, 2}, {1, 0, 0}, {0, 2, 3}, {3, 1, 0},
        {1, 0, 1}, {2, 3, 2}, {2, 3, 3}, {3, 1, 1},
    };

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static V2? CubeDirToFaceDir(int face, V3 dir) {
        V2 ret;
        switch (face>>1) {
            case 0:
                ret = V2.Normalize(new V2(dir.Y, dir.Z));
                break;
            case 1:
                ret = V2.Normalize(new V2(dir.Z, dir.X));
                break;
            case 2:
            default:
                ret = V2.Normalize(new V2(dir.X, dir.Y));
                break;
        }
        return ret.LengthSquared() < 1e-6 ? null : ret;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static V2 CubePosToFacePos(int face, V3 pos) {
        switch (face>>1) {
            case 0:
                return new V2(pos.Y, pos.Z);
            case 1:
                return new V2(pos.Z, pos.X);
            case 2:
            default:
                return new V2(pos.X, pos.Y);
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static V3 FacePosToCubePos(int face, V2 pos) {
        switch(face) {
            case 0:
                return new V3(0.0f, pos.X, pos.Y);
            case 1:
                return new V3(1.0f, pos.X, pos.Y);
            case 2:
                return new V3(pos.Y, 0.0f, pos.X);
            case 3:
                return new V3(pos.Y, 1.0f, pos.X);
            case 4:
                return new V3(pos.X, pos.Y, 0.0f);
            case 5:
            default:
                return new V3(pos.X, pos.Y, 1.0f);
        }
    }
}
