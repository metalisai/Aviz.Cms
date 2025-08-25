namespace Aviz.Cms;

using V3 = System.Numerics.Vector3;

public struct I3 {
    public int x, y, z;

    public I3(int x, int y, int z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
    public V3 ToV3() => new V3(x, y, z);
    public override string ToString() => $"({x} {y} {z})";
    public static I3 operator+(I3 a, I3 b) => new I3(a.x+b.x, a.y+b.y, a.z+b.z);
    public static I3 operator-(I3 a, I3 b) => new I3(a.x-b.x, a.y-b.y, a.z-b.z);
    public static bool operator==(I3 a, I3 b) => a.x == b.x && a.y == b.y && a.z == b.z;
    public static bool operator!=(I3 a, I3 b) => !(a==b);
    public override bool Equals(object? obj) => obj is I3 other && Equals(other);
    public override int GetHashCode() => HashCode.Combine(x, y, z);
};
